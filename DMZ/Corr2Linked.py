import sys

##f_corr
#Gene 1 Gene 2  CorCoeff
#Xelaev18043180m        Xelaev18032367m 0.960043787629994
#Xelaev18043180m        Xelaev18037656m 0.0476190476190476
##f_DB
#Gene stable ID  Gene name
#ENSG00000115204 MPV17   GO:0000002
#ENSG00000198836 OPA1    GO:0000002
##XEN2HS
#XT-MATCH        42sp50.L|gene16313|rna47641     NotAvail        42SP50  42Sp50.L|Xelaev18032367m|100.0
#BOTH    Acod1.L|gene15951|rna46631      ACOD1   ACOD1   irg1.L|Xelaev18013190m|100.0
#BOTH    LOC100036803|gene41709|rna84243 MKNK2   MKNK2   mknk2.S|Xelaev18009863m|96.9

if len(sys.argv) < 3:
    print ("[USAGE] Corr2Linked.py [Corr file] [GeneOntology.clean] [XENLA2HUMAN.DB] ")
    sys.exit(1)

f_corr_name = sys.argv[1]
f_DB_name = sys.argv[2]
f_XEN2HS_name = sys.argv[3]
f_out_name = f_corr_name.split('.txt')[0]+'.90linkage.txt'

try:
    f_corr = open(f_corr_name,'r')
    f_DB = open(f_DB_name,'r')
    f_XEN2HS = open(f_XEN2HS_name,'r')

except:
    print ('File not found!')
    sys.exit(1)

GENE2GO = dict()
f_DB.readline()
for line in f_DB:
    tokens = line.strip().split('\t')
    GeneName = tokens[0]
    GO_name = tokens[1]
    if GeneName not in GENE2GO:
        GENE2GO[GeneName] = set()
    GENE2GO[GeneName].add(GO_name)
f_DB.close()

XEN2HUMAN = dict()
for line in f_XEN2HS:
    tokens = line.strip().split('\t')
    HS_name = tokens[2]
    XEN_name = tokens[4]
    if XEN_name == 'NotAvail' or HS_name == 'NotAvail':
        continue
    XEN_name = XEN_name.split('|')[1]
    if XEN_name not in XEN2HUMAN:
        XEN2HUMAN[XEN_name] = HS_name

f_XEN2HS.close()

f_out = open(f_out_name,'w')
f_corr.readline()
for line in f_corr:
    tokens = line.strip().split('\t')
    GENE_1 = tokens[0]
    GENE_2 = tokens[1]
    Corr_coeff = float(tokens[2])
    if abs(Corr_coeff) >= 0.9:
        if GENE_1 in XEN2HUMAN and GENE_2 in XEN2HUMAN:
            HS_GENE_1 = XEN2HUMAN[GENE_1]
            HS_GENE_2 = XEN2HUMAN[GENE_2]
            if HS_GENE_1 in GENE2GO and HS_GENE_2 in GENE2GO:
                HS_GENE_1_GO = GENE2GO[HS_GENE_1]
                HS_GENE_2_GO = GENE2GO[HS_GENE_2]
                if len(HS_GENE_1_GO & HS_GENE_2_GO) > 0:
                    f_out.write('%s\t%s\t%.20f\tLinked\n'%(GENE_1,GENE_2,Corr_coeff))
                else:
                    f_out.write('%s\t%s\t%.20f\tNot Linked\n'%(GENE_1,GENE_2,Corr_coeff))

f_corr.close()
f_out.close()
