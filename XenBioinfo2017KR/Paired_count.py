
import sys
import re

#BWA
#head of sam file
#@SQ    SN:Kitl|ENSMUST00000105283      LN:5645
#@SQ    SN:Gm9476|ENSMUST00000179109    LN:615

#tail of sam file
#SN840:960:HJTY7BCXY:1:1101:10000:100501        83      Calr|ENSMUST00000003912 197     60      28S73M  =       201     -69     GACTGGAGTTCAGACGTGTGCTCTTCCGATCTATTTCAAAGAGCAGTTCTTGGACGGAGATGCCTGGACCAACCGCTGGGTCGAATCCAAACATAAGTCCG        HIIIIHIIIIIIIIIIHIIIIHIIIIIIIIIIIIIHIIHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD   NM:i:0  MD:Z:73 AS:i:73 XS:i:0

fnm_sam = sys.argv[1]
seq_id2len = dict()
pair_freq = dict()

f_sam = open(fnm_sam,'r')

for line in f_sam:
        if line.startswith('@SQ'):
                tokens = line.strip().split()
                seq_id = re.sub(r'^SN:','',tokens[1])
                seq_len = int(re.sub(r'^LN:','',tokens[2]))
                seq_id2len[seq_id] = seq_len
                continue

        elif line.startswith('@'):
                continue

        tokens = line.strip().split()
        if len(tokens) < 3:
                continue
        pair_id = tokens[2]
        sam_flag = int(tokens[1])
        pair_check = tokens[6]
        frag1_pos = int(tokens[3])
        frag2_pos = int(tokens[7])

        if(pair_id == '*'):
                continue
        if(sam_flag&4):
                continue
        if(frag1_pos == frag2_pos):
                continue

        if(pair_check == '='):
                if not pair_freq.has_key(pair_id):
                        pair_freq[pair_id] = 0
                pair_freq[pair_id] += 1

f_sam.close()
sum_pair_freq = sum(pair_freq.values())

f_out = open(fnm_sam.split('.')[0]+'_RPKM+COUNT.txt','w')
f_out.write('GENE\tCOUNT\tRPKM\n')
for gene_name in sorted(pair_freq.keys()):
        trans_id = gene_name
        trans_read_count = pair_freq[gene_name]
        trans_read_len = seq_id2len[trans_id]
        trans_RPKM = trans_read_count/((trans_read_len/1000.0)*(sum_pair_freq/1000000.0))
        f_out.write('%s\t%d\t%3f\n'%(gene_name,trans_read_count,trans_RPKM))

f_out.close()
