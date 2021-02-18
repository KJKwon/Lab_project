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
tbl.complete  = tbl.complete[,c(-1,-6)]
tbl.male = tbl.complete[tbl.complete$SEX == 'M',]
tbl.female = tbl.complete[tbl.complete$SEX == 'F',]
tmp.select = melt(tbl.complete)
##When considering sex differences
#tmp.select = tmp.select[,c(-1)]
tmp.select$AGE = factor(tmp.select$AGE, levels = sort(unique(tmp.select$AGE)))
tmp.mean = melt(acast(data = tmp.select, AGE~variable, fun = mean))
tbl.mean = tmp.mean[tmp.mean$Var2 == 'PTP4A3',]
tbl.select = tmp.select[tmp.select$variable == 'PTP4A3',]
p = ggplot(tbl.select)+
  geom_boxplot( aes(x = AGE, y = value, group = AGE))+ 
  geom_point(aes(x= AGE, y = value))+
  ggtitle('GTExV8 Ptp4a3')+
  theme(text = element_text(size = 9), legend.title = element_blank(), 
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 9), axis.text.y = element_text(colour = 'black', size = 9),
        axis.title.x = element_text(colour = 'black', size = 9), axis.ticks = element_line(size = 0.9),
        plot.title = element_text(hjust = 0.5, size  = 15)) + 
  ylab('TPM') + xlab('AGE')
  
p
