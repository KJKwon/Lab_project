library(reshape2)
tbl = read.table("TairaUeno_XENLA_tissue.JGIv18pV3_cdna_final.bwa_mem+rsem.indiv_tpm.txt", sep = '\t', header = TRUE, row.names = 1)
threshold = 0.0689519736
tbl.clean = tbl[apply(tbl,1,function(x) median(x) > threshold),]
print('Pearson correlation calculation...')
#Pearson correlation
p.cor = cor(t(tbl.clean), method = 'pearson')
print('Correlation complete!')
p.cor[lower.tri(p.cor, diag = TRUE)] = NA
p.cor = na.omit(melt(p.cor))
print('Duplicate removed!')
colnames(p.cor) = c('Gene1','Gene2','Cor_coeff')
p.cor = p.cor[abs(p.cor$Cor_coeff) > 0.9,]
print('Filtering completed!')
write.table(p.cor, "TairaUeno_XENLA_tissue.JGIv18pV3_cdna_final.bwa_mem+rsem.indiv_tpm_pearson+cor.clean.txt", row.names = FALSE, quote = FALSE, sep = '\t')
#Spearman correlation
print('Spearman correlation calculation...')
s.cor = cor(t(tbl.clean), method = 'spearman')
print('Correlation complete!')
s.cor[lower.tri(s.cor, diag = TRUE)] = NA
s.cor = na.omit(melt(s.cor))
print('Duplicate removed!')
colnames(s.cor) = c('Gene1','Gene2','Cor_coeff')
s.cor = s.cor[abs(s.cor$Cor_coeff) > 0.9,]
print('Filtering completed!')
write.table(s.cor, "TairaUeno_XENLA_tissue.JGIv18pV3_cdna_final.bwa_mem+rsem.indiv_tpm_spearman+cor.clean.txt", row.names = FALSE, quote = FALSE, sep = '\t')
