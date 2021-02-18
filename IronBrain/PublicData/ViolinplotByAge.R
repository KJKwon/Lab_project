library(ggplot2)
library(reshape2)
tbl_tpm = read.table('GTExV8_SN_gene_tpm.gct', header = TRUE, check.names = FALSE)
tbl_tpm = tbl_tpm[,c(-1)]
gene_candi = c('GADD45B','PTP4A3','SULT1A1','VEGFB')
tbl_tpm.selected = tbl_tpm[tbl_tpm$Description %in% gene_candi,]
rownames(tbl_tpm.selected) = tbl_tpm.selected$Description
tbl_tpm.selected = t(tbl_tpm.selected[,c(-1)])
tbl_pheno = read.table('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.full.txt', header = TRUE,row.names = 1)
tbl.complete = merge(tbl_tpm.selected,tbl_pheno, by = intersect(0,0))
tbl.complete  = tbl.complete[,c(-1)]
##When considering sex differences
#tmp.select = tmp.select[,c(-1)]
Candi.gene ='VEGFB'
Select.sex = 'F'
tbl.select = tbl.complete[tbl.complete$SEX == Select.sex,]
#tmp.select = tbl.select
tmp.select = melt(tbl.select)
tmp.select$AGE = factor(tmp.select$AGE, levels = sort(unique(tmp.select$AGE)))
tmp.mean = melt(acast(data = tmp.select, AGE~variable, fun = mean))
tbl.mean = tmp.mean[tmp.mean$Var2 == Candi.gene,]
tbl.select = tmp.select[tmp.select$variable == Candi.gene,]
p = ggplot(tbl.select)+
  geom_violin( aes(x = AGE, y = value, group = AGE))+ 
  geom_jitter(aes(x= AGE, y = value), position = position_jitter(0.1), size = 1)+
  ggtitle(paste('GTExV8','Female','SN',Candi.gene,'expression'))+
  theme(text = element_text(size = 15), legend.title = element_blank(), 
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15), axis.text.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15), axis.ticks = element_line(size = 0.9),
        plot.title = element_text(hjust = 0.5, size  = 20)) + 
  ylab('TPM') + xlab('AGE')

p
