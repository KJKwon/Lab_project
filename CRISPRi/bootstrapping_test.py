import sys
import random
from scipy.stats import spearmanr

def sgRNA2ref(raw_sgRNA_list, sgRNA_DB,sampling_num):
    sgRNA_num2count = dict()
    tmp_sgRNA_sampling = random.choices(raw_sgRNA_list, k = sampling_num)
    for tmp_sgRNA in tmp_sgRNA_sampling:
        if tmp_sgRNA in sgRNA_DB:
            tmp_sgRNA_num = sgRNA_DB[tmp_sgRNA]
            if tmp_sgRNA_num not in sgRNA_num2count:
                sgRNA_num2count[tmp_sgRNA_num] = 0
            sgRNA_num2count[tmp_sgRNA_num] += 1

    return sgRNA_num2count

            
#>@NB552259:10:HFKC7BGXF:1:11101:10001:11472 1:N:0:AATCCAAT
#GGCCCAAGAGGAGAAGCTCG
#>@NB552259:10:HFKC7BGXF:1:11101:10001:8763 1:N:0:AATCCAAT
#CGCGGACAGCGTGGGAAAGA
def fasta2list(file_name):
    f_in = open(file_name,'r')
    tmp_sgRNA_list = []
    for line in f_in:
        if not line.startswith('>'):
            tmp_sgRNA = line.strip()
            tmp_sgRNA_list.append(tmp_sgRNA)
    
    f_in.close()
    return tmp_sgRNA_list

def cor_output(sgRNA1, sgRNA2):
    total_sgRNA_clean = set(sgRNA1.keys()) & set(sgRNA2.keys())
    total_sgRNA_clean = list(total_sgRNA_clean)
    query1_clean = []
    query2_clean = []
    for tmp_sgRNA_num in total_sgRNA_clean:
        query1_clean_sgRNA = sgRNA1[tmp_sgRNA_num]
        query2_clean_sgRNA = sgRNA2[tmp_sgRNA_num]
        query1_clean.append(query1_clean_sgRNA)
        query2_clean.append(query2_clean_sgRNA)
    corr, _ = spearmanr(query1_clean, query2_clean) 
    
    return corr
    
#s_10000000      AAAAAAAAAATACTGAGAGA    GATA3
#s_10000001      AAAAAAAAGAGGAGGGACGG    ANKH
if len(sys.argv) < 2:
    print ('[USAGE] bootstrapping_test.py [FASTA 1] [FASTA 2]')
    sys.exit(1)
else:    
    f_fasta1_nm = sys.argv[1]
    f_fasta2_nm = sys.argv[2]
    f_out_nm = f_fasta1_nm.split('_S')[0] + '_' + f_fasta2_nm.split('_S')[0] +'_correlation_bootstrap_output.txt'
    f_out = open(f_out_nm,'w')
    f_ref = open('broadgpp-dolcetto-targets-set_mageck.txt','r')
    sgRNA2num = dict()
    for line in f_ref:
        token = line.strip().split()
        tmp_sgRNA_seq = token[1]
        tmp_sgRNA_num = token[0]
        sgRNA2num[tmp_sgRNA_seq] = tmp_sgRNA_num

    fasta1_sgRNA = fasta2list(f_fasta1_nm)
    fasta2_sgRNA = fasta2list(f_fasta2_nm)
    f_out.write('Sampling_num\tCorr_coeff\n')
    for i in [0.2,0.25,0.5,1,2,3,4,5]:
        select_num = int(len(fasta1_sgRNA)*i)
        fasta1_ref_sgRNA = sgRNA2ref(fasta1_sgRNA,sgRNA2num,select_num)
        fasta2_ref_sgRNA = sgRNA2ref(fasta2_sgRNA,sgRNA2num,select_num)
        tmp_cor = cor_output(fasta1_ref_sgRNA, fasta2_ref_sgRNA)
        f_out.write('%f\t%.3f\n'%(i,tmp_cor))

f_out.close()
