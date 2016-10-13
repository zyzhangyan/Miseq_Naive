#!/usr/bin/env python
# encoding: utf-8
"""
1.0.py -p project

Created by Mingchen on 2015-05-04.
Copyright (c) 2015 __MyCompanyName__. All rights reserved

"""
import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
import traceback, tempfile
import Bio.Alphabet 
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from misc_prepare_pbs import *
#from bsub import bsub
from mytools import *
try:
    import cPickle as pickle
except ImportError:
    import pickle

### -- BEGIN -- Change fastq file ID
def change_id(prj_folder,ORG_FOLDER,prj_name):
	print 'Changing ID'
 	flag, counts  = 0, 0
	fastq = csv.reader(open(glob.glob('%s/%s/*.fastq'%(prj_folder,ORG_FOLDER))[0]), delimiter = '\t')
	fq = csv.writer(open('%s/%s/%s.fq'%(prj_folder,ORG_FOLDER,prj_name),'w'),delimiter = '\t')
	for index,line in enumerate(fastq):
		if int(index) % 4 == 0:
			flag += 1
			fq.writerow(['@'+str(flag)])
		else:
			fq.writerow(line)
			print flag
### --END -- Change fastq file ID

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def trim_fastq_by_quality_v2(the_file, prj_folder, bad_list):
	handle = open(the_file, "rU")
	fname = retrieve_name_body(the_file)
	print "Triming...",fname
	trim_file = "%s/%s_trimed.fastq"%(prj_tree.origin,fname)
	writer = open(trim_file, "w")
	for record in SeqIO.parse(handle, "fastq") :
		quality_type = list(record.letter_annotations)[0]
		quality_list = record.letter_annotations[quality_type]
		position_list = []
		for index in range(0,len(quality_list)):
			if quality_list[index] > 20:
				position_list.append(index)
		try:
			new_record = record[position_list[0] : position_list[-1]+1]
			if len(new_record.seq) == len(new_record.letter_annotations[quality_type]):
				SeqIO.write(new_record, writer, "fastq")
			else:
				bad_list.append(record.id)
		except:
			bad_list.append(record.id)
			pass
	handle.close()
	writer.close()
	return bad_list

def trim_fastq_by_quality(the_file, prj_folder, bad_list):
	handle = open(the_file, "rU")
	fname = retrieve_name_body(the_file)
	print "Triming...",fname
	trim_file = "%s/%s_trimed.fastq"%(prj_tree.origin,fname)
	writer = open(trim_file, "w")
	for record in SeqIO.parse(handle, "fastq") :
		quality_type = list(record.letter_annotations)[0]
		quality_list = record.letter_annotations[quality_type]
		position_list = []
		for index in range(0,len(quality_list)):
			if quality_list[index] > 20:
				position_list.append(index)
		new_record = record[position_list[0] : position_list[-1]+1]
		SeqIO.write(new_record, writer, "fastq")
	handle.close()
	writer.close()
def unique_fasta(prj_folder):
	handle = "%s/1.2-merged-fastq-file/test0513.extendedFrags.fasta"%prj_folder
	reader = SeqIO.parse(handle, "fasta")
	fname, suffix = os.path.splitext(handle)
	writer = open("%s_unique.fasta"%fname,"w")
	handle_dict, handle_dict_unique, dict_unique = {}, {}, {}
	for index, record in enumerate(reader):
		handle_dict[record.id] = record.seq
		handle_dict_unique.setdefault(record.seq, []).append(record.id)
	for seq, ID in handle_dict_unique.items():
		if len(ID) >= 2:
			print len(ID)
		dict_unique["%s_%d"%(ID[0], len(ID))] = seq
	for ID, seq in dict_unique.items():
		seqrecord = SeqRecord(seq, id =ID)
		SeqIO.write(seqrecord, writer, "fasta")
	print "The number of unique reads in fasta file is %d"%len(dict_unique)

def get_assignment_and_recombanation_info(infile):
	fname = retrieve_name_body(infile)
	count_v,count_d,count_j = 0,0,0 
	c_v, c_d, c_j, result ='', '', '', []
	outfile = open("%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, fname),"w")
	writer = csv.writer(outfile,delimiter = "\t")
	outfile2 = open("%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, fname),"w")
	writer2 = csv.writer(outfile2,delimiter = "\t")
	reader = csv.reader(open(infile,"rU"),delimiter = "\t")
	for line in reader:
		con = str(line)
		con=con.replace('\'','')
		con=con[1:-1]
		hit = re.findall('# Query',con)
		hit_v=re.match(r'^V.+:.+',con)
		hit_d=re.match(r'^D.+:.+',con)
		hit_j=re.match(r'^J.+:.+',con)
		if hit_v:
			a_v=hit_v.group()
			b_v = a_v.split(',')
			if not(c_v == b_v[1]):
				c_v=b_v[1]
				count_v += 1
				writer.writerow(b_v)
		if hit_d:
			a_d=hit_d.group()
			b_d = a_d.split(',')
			if not(c_d == b_d[1]):
				c_d=b_d[1]
				count_d += 1
				writer.writerow(b_d)
		if hit_j:
			a_j=hit_j.group()
			b_j = a_j.split(',')
			if not(c_j == b_j[1]):
				c_j=b_j[1]
				count_j += 1
				writer.writerow(b_j)
		
		if hit:
			con = [':'.join(con.split(":")[1:])]
			con[0].replace(" ", "")
			result.append(con)
		if len(line) >= 7:
			if line[-1] == '+' or line[-1] == '-':
				result[-1].extend(line)
	writer2.writerows(result)

def main():
	print "Begin!"
	#'''
	
	print "Gunzip..."
	try:
		infiles = glob.glob("%s/*.fastq.gz"%(prj_tree.origin))
		infiles = sorted(infiles)

		print infiles[0]
		print infiles[1]	
		os.chdir("%s"%(prj_tree.origin))
		os.system("rm %s/*.fastq"%(prj_tree.origin))
		for infile in infiles:
			gunzip = subprocess.call("gunzip %s "%(infile),shell=True)
	except:
		print "The zip file is not exists!"
		pass
	#'''
	
	#'''
	if os.path.exists("%s/%s.assembled.fastq"%(prj_tree.origin, prj_name)):
		pass
	else:
		print "Merging..."
		infiles = glob.glob("%s/*.fastq"%(prj_tree.origin))
		infiles = sorted(infiles)
		print infiles[0]
		print infiles[1]	
		os.chdir("%s"%(prj_tree.origin))
		merge = subprocess.call("pear -j 4 -q 20 -f %s -r %s -o %s "%(infiles[0],infiles[1], prj_name),shell=True)
	#'''
	
	#'''
	if os.path.exists("%s/%s.assembled_trimed.fastq"%(prj_tree.origin, prj_name)):
		pass
	else:
		print "Quality contorl..."
		infiles = glob.glob("%s/*.assembled.fastq"%(prj_tree.origin))
		trim_files, bad_list = [], []
		for the_file in infiles:
			trim_fastq_by_quality(the_file,prj_folder,bad_list)
	#sys.exit(0)
	#'''	
	'''
	print "Filter..."
	infiles = glob.glob("%s/*assembled_trimed.fastq"%(prj_tree.origin))
	infiles = sorted(infiles)
	r1_infile, r1_id_list = SeqIO.index(infiles[0], "fastq"), []
	r2_infile, r2_id_list = SeqIO.index(infiles[1], "fastq"), []
	
	for ids in r1_infile.values():
		if len(ids.seq) > 10:
			r1_id_list.append(ids.id)
	for	ids in r2_infile.values():
		if len(ids.seq) > 10:
			r2_id_list.append(ids.id)
	pair_reads = set(r1_id_list) & set(r2_id_list)
	
	r1_file_writer, r2_file_writer = open("%s/%s_R1_filter.fastq"%(prj_tree.origin, prj_name), 'w'), open("%s/%s_R2_filter.fastq"%(prj_tree.origin, prj_name), 'w')
	for read_id in list(pair_reads):
		if read_id not in bad_list:
			SeqIO.write(r1_infile[read_id], r1_file_writer, "fastq")
			SeqIO.write(r2_infile[read_id], r2_file_writer, "fastq")
	'''
	
	#'''
	if os.path.exists("%s/%s.assembled_trimed.fasta"%(prj_tree.origin, prj_name)):
		pass
	else:
		print "Convert fastq to fasta..."
		merged_file = "%s/%s.assembled_trimed.fastq"%(prj_tree.origin, prj_name)
		fname, suffix = os.path.splitext(merged_file)
		count = SeqIO.convert(merged_file, "fastq","%s.fasta"%fname, "fasta")
		print count
		print "There are  %i records have been Converted!" %(count)
	#'''
	'''
	print "Unique fasta file..."
	unique_fasta(prj_folder)
	
	'''
	#Step 2: Split to little files
	print "Step 2: Split to little files"
	record_iter = SeqIO.parse(open("%s/%s.assembled_trimed.fasta"%(prj_tree.origin, prj_name)), "fasta")
	for i, batch in enumerate(batch_iterator(record_iter, 10000)) :
		filename = "%s/%s_%i.fasta" % (prj_tree.split,prj_name, i+1)
		handle = open(filename, "w")
		count = SeqIO.write(batch, handle, "fasta")
		handle.close()
		print "Wrote %i records to %s" % (count, filename)
	files_num = i+1
	
	#'''
	#Step 5: Mapping, Multiple processing
	print "Begin IgBLAST..."
	cmd = "cp -r /zzh_gpfs02/yanmingchen/HJT-PGM/PGM_UMI_121_20160624/Igblast_database/* %s/"%(prj_tree.igblast_database)
	gunzip = subprocess.call(cmd,shell=True)
	prepare_IgBLAST_jobs(prj_name, prj_tree)
	IgBLAST_jobs = glob.glob("%s/IgBLAST_*.sh" %(prj_tree.jobs))
	IgBLAST_pool = Pool()
	
	#'''#Cluster PBS
	IgBLAST_jobs_ids = IgBLAST_pool.map_async(bsub_jobs, IgBLAST_jobs).get(120)
	print "IgBLAST_jobs: %s IgBLAST_job has been submited."%len(IgBLAST_jobs_ids)
	check_jobs_done(prj_name, prj_tree, "IgBLAST", IgBLAST_jobs_ids)
	print "Waiting for all subprocesses done..."
	IgBLAST_pool.close()
	IgBLAST_pool.join()
	print 'All subprocesses done.'
	#'''
	
	'''#One machine
	IgBLAST_jobs_ids = IgBLAST_pool.map_async(processing_jobs, IgBLAST_jobs).get(120)
	print "IgBLAST_jobs: %s IgBLAST_job has been submited."%len(IgBLAST_jobs_ids)
	print "Waiting for all subprocesses done..."
	IgBLAST_pool.close()
	IgBLAST_pool.join()
	print 'All subprocesses done.'
	'''
	#'''
	
	#"""
	igblast_result_files = glob.glob("%s/IgBLAST_result_*.txt"%(prj_tree.igblast_data))
	pool = Pool()
	for igblast_result_file in igblast_result_files:
		#get_assignment_and_recombanation_info(igblast_result_file)
		pool.apply_async(get_assignment_and_recombanation_info, args=(igblast_result_file,))
	print "Waiting for all subprocesses done..."
	pool.close()
	pool.join()
	print 'All subprocesses done.'
	os.system("cat %s/IgBLAST_result_*_get_assignment_info.txt > %s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_tree.igblast_data, prj_name))
	os.system("cat %s/IgBLAST_result_*_get_recombanation_info.txt > %s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_tree.igblast_data, prj_name))
	
	#"""

if __name__ == '__main__':
	print 'Parent process %s'%os.getpid()
	prj_folder = os.getcwd()
	prj_tree = create_folders(prj_folder)
	prj_tree = ProjectFolders(prj_folder)
	
	prj_name = fullpath2last_folder(prj_tree.home)
	pool_size = multiprocessing.cpu_count()
	main()
