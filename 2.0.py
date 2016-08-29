#!/usr/bin/env python
# encoding: utf-8
"""
2.0-process-igblast-output.py 

Created by Mingchen on 2015-05-15.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

import glob, os, sys, subprocess, re
from Bio import SeqIO
from mytools_ymc import *

def process_VDJ_alignment(choosed_infile,IGBLAST_output_files, germline_type):
	reader_choosed = csv.reader(open(choosed_infile,"rU"),delimiter = "\t")
	choosed_list = []
	for line in reader_choosed:
		choosed_list.append(line.split('*')[0])
	NoMut_alignment_dict, alignment_dict, subject_result_list, subject_result_family_list, NoMut_result_list = {}, {}, [], [], []
	for infile in IGBLAST_output_files:
		print "Processing %s " %infile
		reader = csv.reader(open(infile,"rU"),delimiter = "\t")
		for line in reader:
			if len(line) == 30 and float(line[3]) == float(100.00) and float(line[12]) < 1e-15 and line[0] == germline_type and int(line[10]) == 1 and line[11] == line[15]:
				query_id = line[1]
				NoMut_result_list.append(line)
				NoMut_alignment_dict.setdefault(line[1],[]).append([line[2],line[13]])
			if len(line) == 30 and line[3] >= 80 and float(line[12]) < 1e-15 and line[0] == germline_type and line[2] in choosed_list:
				query_id = line[1]
				alignment_dict.setdefault(line[1],[]).append([line[2],line[13]])
	subject_result_set = set(subject_result_list)
	return NoMut_alignment_dict, alignment_dict, NoMut_result_list


def def_score_function(alignment_dict):
	subject_id_value_dict, value_list = {}, []
	for value in alignment_dict.values():
		#print value
		for item in value:
			#processed_item = (item[0],float(1)/len(value))
			#processed_item = (item[0],1)#simple way
			processed_item = (item[0],float(item[1])/len(value))
			#processed_item = (item[0],float(item[1]))
			value_list.append(processed_item)
	return value_list


def calculate_score(value_list):
	unique_ID = []
	for item in value_list:
		flag=1
		for i in range(0,len(unique_ID)):
			if item[0] == unique_ID[i][0]:
				unique_ID[i][1] += item[1]
				flag=0
				break
		if flag==1:
			unique_ID.append(list(item))
	return unique_ID

def main():
	print "Begin!"
	prj_folder = os.getcwd()
	choosed_infile = "%s/choosed_gene.txt"%(prj_folder)
	IGBLAST_output_files = glob.glob("%s/1.4-IgBLAST-output/IgBLAST_result_*"%(prj_folder))
	
	os.chdir("%s/2.1-germline-score/"%(prj_folder))
	for germline_type in ('V','D','J'):
		NoMut_alignment_dict, alignment_dict, NoMut_result_list = process_VDJ_alignment(choosed_infile,IGBLAST_output_files, germline_type)
		value_list = def_score_function(alignment_dict)
		
		NoMut_line_writer = csv.writer(open("NoMut_line_germeline_%s_score_bitscore.txt"%germline_type,"wt"), delimiter = "\t")
		NoMut_line_writer.writerows(NoMut_result_list)
		
		unique_ID = calculate_score(value_list)
		writer = csv.writer(open("germeline_%s_score_bitscore.txt"%germline_type,"wt"), delimiter = "\t")
		
		NoMut_value_list = def_score_function(NoMut_alignment_dict)
		NoMut_unique_ID = calculate_score(NoMut_value_list)
		NoMut_writer = csv.writer(open("NoMut_germeline_%s_score_bitscore.txt"%germline_type,"wt"), delimiter = "\t")
		
		for item in sorted(unique_ID):
			writer.writerow(item)
		
		for item in sorted(NoMut_unique_ID):
			NoMut_writer.writerow(item)

if __name__ == '__main__':

	create_folders(os.getcwd())
	
	main()
	print "Finished"
