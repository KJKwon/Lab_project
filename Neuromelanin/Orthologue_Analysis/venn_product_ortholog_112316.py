import sys
F_A = open('Human_sig_gene_output.txt','r')
F_B = open('Mouse_sig_output_111716.txt','r')

h_pos_dat = dict()
h_neg_dat = dict()
m_pos_dat = dict()
m_neg_dat = dict()

def Test_out(query,db_comp,db_orth,DB,code):
	identify = code
	scatter_out = dict()
	together_gene = []
	orth_gene = []
	only_gene = []
	not_found_gene = []
	num = 0
	for query_gene in query:
		query_gene_pval = query[query_gene][1]
		query_gene_FC =query[query_gene][0]
		if query_gene_pval < 0.05:
			num += 1
			if query_gene in db_orth:
				orth_gene.append(query_gene)
				if len(db_orth[query_gene]) == 1:
					db_candidate = db_orth[query_gene][0]
					if db_candidate in DB:
		                               	if db_candidate in db_comp:
							scatter_out[query_gene] = [query[query_gene][0],db_comp[db_candidate][0]]
							if db_comp[db_candidate][1] < 0.05 and query[query_gene][0] * db_comp[db_candidate][0] > 0:
								together_gene.append(query_gene)
							else:
								only_gene.append(query_gene)
						else:
							only_gene.append(query_gene)
					else:
						not_found_gene.append(query_gene)
				else:
					check = []
					for i in range(len(db_orth[query_gene])):
						db_candidate = db_orth[query_gene][i]
						if db_candidate in DB:
							if db_candidate in db_comp:
								if db_comp[db_candidate][1]<0.05 and query[query_gene][0] * db_comp[db_candidate][0] > 0:
									scatter_out[query_gene] = [query[query_gene][0],db_comp[db_candidate][0]]
									together_gene.append(query_gene)
									check.append(2)
									break
								else:
									temp_scatter = db_candidate
									check.append(1)
							else:
								check.append(1)
						else:
							continue
					if 2 not in check:
						if 1 in check:
							only_gene.append(query_gene)
						else:
							not_found_gene.append(query_gene)
			else:
				not_found_gene.append(query_gene)
		final_only_gene = []
		if identify[1] == 'A':
			for sig_gene in only_gene:
				if query[sig_gene][0] > 2:
					final_only_gene.append(sig_gene)
		else:
			for sig_gene in only_gene:
				if query[sig_gene][0] < -2:
					final_only_gene.append(sig_gene)
	print num
	
	return together_gene,final_only_gene,orth_gene,identify

def out_write(data_list):
	together = data_list[0]
	only = data_list[1]
	orth = data_list[2]
	name_check = data_list[3]
	name = ''
	if name_check[0] == '1':
		name += 'human'
		if name_check[1] == 'A':
			name+= '_positive_'
		else:
			name+= '_negative_'
	else:
		name += 'mouse'
		if name_check[1] == 'A':
			name+= '_positive_'
		else:
			name+= '_negative_'
	write_list = [together,only,orth]
	for i in range(len(write_list)):
		num = 0
		DB = write_list[i]
		if i == 0:
			f_out = open(name+'together.txt','w')
		elif i == 1:
			f_out = open(name+'only.txt','w')
		else:
			f_out = open(name+'orth.txt','w')
		for gene in DB:
			f_out.write('%s\n'%(gene))
			num += 1
		f_out.write('%d'%(num))
		f_out.close()
	return set(together)
			
F_A.readline()
human_DB = []
for line in F_A:
	token = line.strip().split('\t')
	h_gene = token[0].split('"')[1]
	h_FC = float(token[1])
	h_pval = float(token[5])
	human_DB.append(h_gene.upper())
	if h_FC > 0:
		h_pos_dat[h_gene.upper()] = [h_FC,h_pval]
	else:
		h_neg_dat[h_gene.upper()] = [h_FC,h_pval]
F_A.close()
F_B.readline()
mouse_DB = []
for line in F_B:
	token = line.strip().split('\t')
	m_gene = token[0]
	m_FC = float(token[1])
	m_pval = float(token[2])
	mouse_DB.append(m_gene.upper())
	if m_FC > 0:
		m_pos_dat[m_gene.upper()] = [m_FC,m_pval]
	else:
		m_neg_dat[m_gene.upper()] = [m_FC,m_pval]
F_B.close()

F_orth = open('human_mouse_orth_data.txt','r')
F_orth.readline()
h_m_transfer = dict()
m_h_transfer = dict()
for line in F_orth:
	token = line.strip().split()
	h_orth = token[0]
	m_orth = token[1]
	if h_orth not in h_m_transfer:
		h_m_transfer[token[0].upper()] = [token[1].upper()]
	else:
		h_m_transfer[token[0].upper()].append(token[1].upper())
	if m_orth not in m_h_transfer:
		m_h_transfer[token[1].upper()] = [token[0].upper()]
	else:
		m_h_transfer[token[1].upper()].append(token[1].upper())

F_orth.close()
#human = 1, mouse = 2, positive =A, negative = B
hp_mp = out_write(Test_out(h_pos_dat,m_pos_dat,h_m_transfer,mouse_DB,'1A'))
hn_mn = out_write(Test_out(h_neg_dat,m_neg_dat,h_m_transfer,mouse_DB,'1B'))
mp_hp = out_write(Test_out(m_pos_dat,h_pos_dat,m_h_transfer,human_DB,'2A'))
mn_hn = out_write(Test_out(m_neg_dat,h_neg_dat,m_h_transfer,human_DB,'2B'))

f_pos = open('total_positive_list','w')
f_neg = open('total_negative_list','w')

p_total = hp_mp | mp_hp
n_total = hn_mn | mn_hn
p_total = list(p_total)
n_total = list(n_total)
for g_pos in p_total:
	f_pos.write('%s\n'%(g_pos))
f_pos.write('%d\n'%(len(p_total)))
for g_neg in n_total:
	f_neg.write('%s\n'%(g_neg))
f_neg.write('%d\n'%(len(n_total)))

f_pos.close()
f_neg.close()
print len(h_pos_dat)
print len(h_neg_dat)
print len(m_pos_dat)
print len(m_neg_dat)
