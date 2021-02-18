import sys
###GTExV8_SN_SampleIDList.txt
##SAMPID
#GTEX-11DZ1-0011-R2a-SM-DNZZM
#GTEX-11EMC-0011-R2b-SM-DO114

##GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct
#Name    Description     GTEX-1117F-0226-SM-5GZZ7   ... GTEX-1117F-0226-SM-5GZZ9
#ENSG00000227232.5       WASH7P  187    ...     109

f_list = open('GTExV8_SN_SampleIDList.txt','r')
f_list.readline()
SampleList = []
for line in f_list:
    tmp_ID = line.strip()
    SampleList.append(tmp_ID)

f_list.close()

f_count_tbl = open('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct','r')
f_count_tbl.readline()
f_count_tbl.readline()
Samples = f_count_tbl.readline().strip().split()[2:]
TrueSamples = []
SampleIndex = []
for tmp_TargetSample in SampleList:
    if tmp_TargetSample in Samples:
        tmp_Index = Samples.index(tmp_TargetSample)
        SampleIndex.append(tmp_Index)
        TrueSamples.append(tmp_TargetSample)
Gene2Count = dict()
for line in f_count_tbl:
    tokens = line.strip().split('\t')
    GeneName = tokens[0]
    GeneSymbol = tokens[1]
    CountList = tokens[2:]
    Gene2Count[(GeneName,GeneSymbol)] = [CountList[x] for x in SampleIndex]

f_count_tbl.close()
f_out = open('GTExV8_SN_gene_tpm.gct','w')
f_out.write('Name\tDescription\t%s\n'%('\t'.join(TrueSamples)))
for tmp_Gene in Gene2Count.keys():
    tmp_Count = Gene2Count[tmp_Gene]
    tmp_GeneName = tmp_Gene[0]
    tmp_GeneSymbol = tmp_Gene[1]
    f_out.write('%s\t%s\t%s\n'%(tmp_GeneName,tmp_GeneSymbol,'\t'.join(tmp_Count)))

f_out.close()
