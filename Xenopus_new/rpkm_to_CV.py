import numpy as np
import sys
import os
import math

if len(sys.argv) < 2:
	print 'Usage:[rpkm_to_CV.py][file_name]'
	sys.exit(1)
file_temp = sys.argv[1]
for file_nm in os.listdir('.'):
	if file_temp == file_nm:
		file_name = file_temp
gene_to_CV = dict()
f_in = open(file_name,'r')
f_in.readline()
for line in f_in:
	tokens = line.strip().split('\t')
	gene_name = tokens[0]
	B_list = tokens[2:12]
	C_list = tokens[14:24]
	BC_mean_list = []
	for num in range(len(B_list)):
		BC_mean = (float(B_list[num])+float(C_list[num]))/2
		BC_mean_list.append(BC_mean)
	CV_mean = np.mean(BC_mean_list)
	CV_std = np.std(BC_mean_list)
	CV = CV_std/CV_mean
	gene_to_CV[gene_name] = CV
f_in.close()

f_out = open(file_name.split('.')[0]+'_CV.txt','w')
for gene in gene_to_CV:
	f_out.write('%s\t%f\n'%(gene,math.log(gene_to_CV[gene])))
f_out.close()
