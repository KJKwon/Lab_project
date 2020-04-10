import sys
import numpy as np

#Corcoeff cut off 0.8
if len(sys.argv) < 2:
    print ("[USAGE] tbl2corcoeff.py [data_table]")
    sys.exit(1)

f_in_nm = sys.argv[1]
f_in = open(f_in_nm,'r')

f_in.readline()
Total_ID = []
Total_exp = []
for line in f_in:
    tokens = line.strip().split('\t')
    XenID = tokens[0].split('|')[1]
    ExprVal = [float(x) for x in tokens[1:]]
    temp_exp = np.array(ExprVal)
#First quarter of total expression
    if np.mean(temp_exp) > 0.036:
        Total_ID.append(XenID)
        Total_exp.append(ExprVal)

f_in.close()

f_out_nm = f_in_nm.split('.txt')[0] + '_cor.clean.txt'
f_out = open(f_out_nm,'w')
f_out.write('Gene1\tGene2\tCorrcoeff\n')
for i in range(len(Total_ID)-1):
    for j in range(i+1,len(Total_ID)):
        Gene1 = Total_ID[i]
        Gene2 = Total_ID[j]
        Gene1_exp = Total_exp[i]
        Gene2_exp = Total_exp[j]
        corcoeff = np.corrcoef(Gene1_exp,Gene2_exp)[0][1]
        if abs(corcoeff) >= 0.8:
            f_out.write('%s\t%s\t%.20f\n'%(Gene1,Gene2,corcoeff))

f_out.close()
