import sys
##Fastq
#@FS10000436:43:BPL20321-2819:1:1101:10180:1000 1:N:0:TTGAACCG
#GCTTGTGGAAAGGACGAAACACCGGCCCAGGGGCTGCCGTCCCGGTTTCAGAGCTACAGCAGAAATGCTGTAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTCAAGCTTGGC
#+
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#Adapter = 'TTGTGGAAAGGACGAAACACCG'
if len(sys.argv) < 2:
    print('[USAGE] fastq2stag+sgRNA.py [fastq]')
    sys.exit(1)
else:
    f_in_nm = sys.argv[1]
    f_in = open(f_in_nm,'r')

sgRNA2stag = dict()
ID2sgRNA = dict()
PASS = False
staggers = ['GAAGACCC','ACGCAAC','TGCACC','CAAC','AGC','GC','C']
stag_count = [0,0,0,0,0,0,0,0]
adapter = 'TTGTGGAAAGGACGAAACACCG'
f_adapter = 'TTGTGGAAAGG'
b_adapter = 'ACGAAACACCG'
for line in f_in:
    if line.startswith('@'):
        ID = line.strip()
        PASS = True
    else:
        if PASS == True:
            sequence = line.strip()
            if adapter in sequence:
                adapt_start = sequence.index(adapter)
            elif f_adapter in sequence:
                adapt_start = sequence.index(f_adapter)
            elif b_adapter in sequence:
                adapt_start = sequence.index(b_adapter)-10
            else:
                continue
            sgRNA_start = adapt_start + 22
            sgRNA = sequence[sgRNA_start:sgRNA_start+20]
            if adapt_start != 0:
                stagger_region = sequence[:adapt_start]
                for i in range(len(staggers)):
                   tmp_stagger = staggers[i]
                   if len(tmp_stagger) <= len(stagger_region):
                        if tmp_stagger == stagger_region[(-len(tmp_stagger)):]:
                            stag_count[i] += 1
                            final_stag = tmp_stagger
                            if sgRNA not in sgRNA2stag:
                                sgRNA2stag[sgRNA] = [0,0,0,0,0,0,0,0]
                            sgRNA2stag[sgRNA][i] += 1
                            ID2sgRNA[ID] = sgRNA
                            break
            else:
                ID2sgRNA[ID] = sgRNA
                final_stag = ''
                stag_count[-1] += 1
                if sgRNA not in sgRNA2stag:
                    sgRNA2stag[sgRNA] = [0,0,0,0,0,0,0,0]
                sgRNA2stag[sgRNA][-1] += 1

            PASS = False

ID_list = sorted(list(ID2sgRNA.keys()))
f_stag_out = open(f_in_nm.split('.txt')[0] + '_stag+sgRNA.txt','w')
f_stag_out.write('##Total_stag_out = %d\t'%(sum(stag_count))+'GAAGACCC = {}\tACGCAAC = {}\tTGCACC = {}\tCAAC = {}\tAGC = {}\tGC = {}\tC = {}\tNA = {}\n'.format(*stag_count))
f_stag_out.write('sgRNA\tstag_total\t%s\n'%('\t'.join(staggers)))
f_sgRNA_out = open(f_in_nm.split('.txt')[0] + '_sgRNA.fasta','w')
for tmp_ID in ID_list:
    tmp_sgRNA = ID2sgRNA[tmp_ID]
    f_sgRNA_out.write('>%s\n%s\n'%(tmp_ID,tmp_sgRNA))

f_sgRNA_out.close()

sgRNA_list = sorted(list(sgRNA2stag.keys()))
for tmp_sgRNA in sgRNA_list:
    tmp_stag_list = sgRNA2stag[tmp_sgRNA]
    sum_stag = sum(tmp_stag_list)
    tmp_stag_list = [str(x) for x in tmp_stag_list]
    f_stag_out.write('%s\t%d\t%s\n'%(tmp_sgRNA,sum_stag,'\t'.join(tmp_stag_list)))

f_stag_out.close()
                 
