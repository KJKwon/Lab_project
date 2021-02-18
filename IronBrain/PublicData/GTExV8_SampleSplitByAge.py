import sys
##GTExV8_SN_gene_reads.gct
#Name    Description     GTEX-11DZ1-0011-R2a-SM-DNZZM ... GTEX-11EMC-0011-R2b-SM-DO114
#ENSG00000223972.5       DDX11L1 0 ... 0

##GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
#SUBJID  SEX     AGE     DTHHRDY
#GTEX-1117F      2       60-69   4
#GTEX-111CU      1       50-59   0
#SUBJID is specific for each specific region

def WriteFile(FileName,IDtoCOUNT,GeneTotalList,IDList):
    f_out = open(FileName,'w')
    f_out.write('Name\tDescription\t%s\n'%('\t'.join(IDList)))
    for tmp_GeneName in GeneTotalList:
        tmp_GeneCount = []
        for tmp_ID in IDList:
            tmp_Count = ID2COUNT[tmp_ID][tmp_GeneName]
            tmp_GeneCount.append(tmp_Count)
        f_out.write('%s\t%s\n'%('\t'.join(list(tmp_GeneName)),'\t'.join(tmp_GeneCount)))
    f_out.close()

f_count = open('GTExV8_SN_gene_tpm.gct','r')
HEAD2ID = dict()
ID2COUNT = dict()
GeneNameList = []
SampleList = f_count.readline().strip().split()[2:]
for tmp_sample in SampleList:
    tmp_SampleHead = 'GTEX-'+tmp_sample.split('-')[1]
    if tmp_SampleHead not in HEAD2ID:
        HEAD2ID[tmp_SampleHead] = []
    HEAD2ID[tmp_SampleHead].append(tmp_sample)

for line in f_count:
    tokens = line.strip().split('\t')
    tmp_GeneName = tuple(tokens[:2])
    GeneNameList.append(tmp_GeneName)
    tmp_Count = tokens[2:]
    for i in range(len(SampleList)):
        tmp_Sample = SampleList[i]
        if tmp_Sample not in ID2COUNT:
            ID2COUNT[tmp_Sample] = dict()
        ID2COUNT[tmp_Sample][tmp_GeneName] = tmp_Count[i]

f_count.close()

f_pheno = open('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt','r')
PHENO2ID = dict()
f_pheno.readline()
for line in f_pheno:
    tokens = line.strip().split('\t')
    tmp_IDHead = tokens[0]
    tmp_sex = int(tokens[1])
    tmp_age = tokens[2]
    if tmp_sex == 1:
        tmp_sex = 'male'
    else:
        tmp_sex = 'female'
    if tmp_IDHead in HEAD2ID:
        if tmp_age not in PHENO2ID:
            PHENO2ID[tmp_age] = dict()
        if tmp_sex not in PHENO2ID[tmp_age]:
            PHENO2ID[tmp_age][tmp_sex] = []
        PHENO2ID[tmp_age][tmp_sex] += HEAD2ID[tmp_IDHead]

f_pheno.close()

for tmp_age in PHENO2ID.keys():
    tmp_sex = 'male'
    tmp_IDList = PHENO2ID[tmp_age][tmp_sex]
    tmp_FileName = 'GTEx_Analysis_v8_SN_GeneCount_'+tmp_age+'_'+tmp_sex+'.txt'
    WriteFile(tmp_FileName, ID2COUNT,GeneNameList,tmp_IDList)
    tmp_sex = 'female'
    tmp_IDList = PHENO2ID[tmp_age][tmp_sex]
    tmp_FileName = 'GTEx_Analysis_v8_SN_GeneCount_'+tmp_age+'_'+tmp_sex+'.txt'
    WriteFile(tmp_FileName, ID2COUNT,GeneNameList,tmp_IDList)
