#Transcript stable ID    Gene stable ID
#ENSRNOT00000041720      ENSRNOG00000031780
##Count File
#geneID Rat_15Mo_3L     Rat_15Mo_3R     Rat_15Mo_4L     Rat_15Mo_4R     Rat_6Mo_5L      Rat_6Mo_5R      Rat_6Mo_6L      Rat_6Mo_6R
#1700020D05Rik|ENSRNOT00000084065       6       4       6       17      13      3       7       15
if len(sys.argv) < 3:
    print ("[USAGE] TranscriptID2GeneID.py [TranscriptID2GeneID ensembl] [Count data]")
    sys.exit(1)
else:
    f_DB_name = sys.argv[1]
    f_count_name = sys.argv[2]
    if f_DB_name not in os.listdir('.') or f_count_name not in os.listdir('.'):
        print ('File not found!')
        sys.exit(1)
    else:
        f_DB = open(f_DB_name,'r')
        f_count = open(f_count_name,'r')

f_DB.readline()
Trans2Gene = dict()
for line in f_DB:
    tokens = line.strip().split('\t')
    TransID = tokens[0]
    GeneID = tokens[1]
    Trans2Gene[TransID] = GeneID

f_DB.close()

f_out = open(f_count_name.split('.')[0]+'_GeneID.txt','w')
f_out.write(f_count.readline())
GeneID2Count = dict()
for line in f_count:
    tokens = line.strip().split('\t')
    TransID = tokens[0].split('|')[1]
    Counts = tokens[1:]
    GeneID = Trans2Gene[TransID]
    if GeneID not in GeneID2Count:
        GeneID2Count[GeneID] = Counts
    else:
        new_count = []
        old_count = GeneID2Count[GeneID]
        for i in range(len(Counts)):
            tmp_count = int(old_count[i]) + int(Counts[i])
            new_count.append(str(tmp_count))
        GeneID2Count[GeneID] = new_count

f_count.close()

GeneID_list = sorted(GeneID2Count.keys())
for GeneID in GeneID_list:
    f_out.write('%s\t%s\n'%(GeneID,'\t'.join(GeneID2Count[GeneID])))

f_out.close()
