import sys
import os
#logFC   logCPM  PValue  FDR
#ENSRNOG00000060048      12.6683572494449        5.23937555192889        1.56401184633849e-154   2.09374265869333e-150

if len(sys.argv) < 3:
    print ("[USAGE] edgeR_GeneID2GeneName.py [GeneID2GeneName ensembl] [edgeR output]")
    sys.exit(1)
else:
    f_DB_name = sys.argv[1]
    f_edgeR_name = sys.argv[2]
    if f_DB_name not in os.listdir('.'):
        print ('File not found!')
        sys.exit(1)
    else:
        f_DB = open(f_DB_name,'r')
        f_edgeR_out = open(f_edgeR_name,'r')

geneID2gene = dict()
f_DB.readline()
for line in f_DB:
    token = line.strip().split('\t')
    gene_nm = token[1]
    geneID = token[0]
    geneID2gene[geneID] = gene_nm
    
f_DB.close()
f_total_out = open('15Mo34_vs_6Mo56.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.edgeR.GeneSymbol.clean+cpm.txt','w')
#f_clean_out = open('15Mo34_vs_6Mo56.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.edgeR.GeneSymbol.edgeR.clean.txt','w')
edgeR_ttl = f_edgeR_out.readline()
f_total_out.write('GeneName\t'+ edgeR_ttl)
#f_clean_out.write('GeneName\t'+ edgeR_ttl)

geneID_list = list(geneID2gene.keys())
for line in f_edgeR_out:
    tokens = line.strip().split('\t')
    geneInfo = tokens[1:]
    geneID = tokens[0]
    if geneID in geneID_list:
        gene_name = geneID2gene[geneID]
    else:
        gene_name = 'NA'
    f_total_out.write('%s\t%s\n'%(gene_name,'\t'.join(geneInfo)))
#    logFC = float(tokens[1])
#    FDR = float(tokens[-1])
#    if abs(logFC) > 1 and FDR < 0.05:
#        f_clean_out.write('%s\t%s\n'%(gene_name,'\t'.join(geneInfo)))

f_edgeR_out.close()
f_total_out.close()
#f_clean_out.close()
