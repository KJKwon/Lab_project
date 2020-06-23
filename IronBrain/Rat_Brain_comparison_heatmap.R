tbl = read.table('SH-SY5Y_IronTreat_Ctrlvs10mM.edgeRQLF_output_clean_RatOrtho.txt', sep = '\t', header = TRUE)
Rat_ortho = unique(tbl$RatOrtho)
#Rat_ortho_up = tbl$RatGene[tbl$Trend == 'IronRich Up']
#Rat_ortho_down = tbl$RatGene[tbl$Trend == 'IronRich down']
Rat_count_tbl = read.table('Rat_Brain_15M_6M_6W_count_Left+Right_GeneName.txt', sep = '\t', header = TRUE, row.names = 1)
Rat_rpkm_tbl = read.table('Rat_Brain_15M_6M_6W_rpkm_Left+Right_GeneName.txt', sep = '\t', header = TRUE , row.names = 1)
Rat_rpkm_candi = Rat_rpkm_tbl[rownames(Rat_rpkm_tbl) %in% Rat_ortho,]
Rat_rpkm_male = Rat_rpkm_candi[,1:8]
Rat_rpkm_female = Rat_rpkm_candi[,9:18]
#Filtering and correlation check
select = apply(Rat_rpkm_female, 1, function(x) sum(x >= 1) >= 5 )
Rat_rpkm_female.clean = Rat_rpkm_female[select,]
tmp_cor = dist(1 - cor(as.matrix(t(Rat_rpkm_female.clean.logFC)), method = 'spearman'))
tmp_clust = hclust(tmp_cor, method = 'average')
tmp_clust = as.dendrogram(tmp_clust)
sample_cor = dist(1 - cor(as.matrix(Rat_rpkm_female.clean.logFC), method = 'spearman'))
sample_clust = hclust(sample_cor, method = 'average')
sample_clust = as.dendrogram(sample_clust)

library(gplots)
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)
#RPKM
tmp_logFC = as.vector(apply(Rat_rpkm_female.clean[,6:10], 1, mean))
Rat_rpkm_female.clean.logFC = Rat_rpkm_female.clean[,1:5] / tmp_logFC
select = apply(Rat_rpkm_female.clean.logFC, 1, function(x) sum( x == 0 ) == 0 )
Rat_rpkm_female.clean.logFC = Rat_rpkm_female.clean.logFC[select,]
Rat_rpkm_female.clean.logFC = Rat_rpkm_female.clean.logFC[apply(Rat_rpkm_female.clean.logFC, 1, function(x) sum(abs(log(x,base = 2)) < 0.3) == 0),]
Rat_rpkm_female.clean.logFC = Rat_rpkm_female.clean.logFC[order(apply(Rat_rpkm_female.clean.logFC[,c(1,2,5)], 1, mean)),]
k = heatmap.2(as.matrix(log(Rat_rpkm_female.clean.logFC, base = 2)), trace = "none", col = Colors,Rowv = FALSE, Colv = FALSE, main = "Female 6Mo vs 6Wo heatmap [Age]")

Rat_rpkm_female.revise = Rat_rpkm_female.clean[rownames(Rat_rpkm_female.clean) %in% rownames(Rat_rpkm_female.clean.logFC),]
Rat_tpm_female.revise = Rat_tpm_female.revise[order(apply(Rat_tpm_female.revise[,1:5], 1, mean)),]
select = apply(Rat_rpkm_female.revise, 1, function(x) sum( x == 0 ) == 0 )
Rat_rpkm_female.revise = Rat_rpkm_female.revise[select,]
Rat_tpm_female.revise = t(t(Rat_rpkm_female.revise) / as.vector(apply(Rat_rpkm_female.revise, 2, sum))) * (10^6)
k = heatmap.2(as.matrix(log(Rat_tpm_female.revise, base = 2)), trace = "none", col = Colors,Rowv = FALSE, Colv = FALSE, main = "Female 6Mo vs 6Wo heatmap [Age]")

##TPM
Rat_tpm_female.clean = t(t(Rat_rpkm_female.clean) / as.vector(apply(Rat_rpkm_female.clean, 2, sum))) * (10^6)
tmp_logFC = as.vector(apply(Rat_tpm_female.clean[,6:10], 1, mean))
Rat_tpm_female.clean.logFC = Rat_tpm_female.clean[,1:5] / tmp_logFC
select = apply(Rat_tpm_female.clean.logFC, 1, function(x) sum( x == 0 ) == 0 )
Rat_tpm_female.clean.logFC = Rat_tpm_female.clean.logFC[select,]
Rat_tpm_female.clean.logFC = Rat_tpm_female.clean.logFC[apply(Rat_tpm_female.clean.logFC, 1, function(x) sum(abs(log(x,base = 2)) < 0.2) == 0),]
Rat_tpm_female.clean.logFC = Rat_tpm_female.clean.logFC[order(apply(Rat_tpm_female.clean.logFC[,1:5], 1, mean)),]
k = heatmap.2(as.matrix(log(Rat_tpm_female.clean.logFC, base = 2)), trace = "none", col = Colors,Rowv = FALSE, Colv = FALSE, main = "Female 6Mo vs 6Wo heatmap [Age]")

Rat_tpm_female.revise = Rat_tpm_female.clean[rownames(Rat_tpm_female.clean) %in% rownames(Rat_tpm_female.clean.logFC),]
Rat_tpm_female.revise = Rat_tpm_female.revise[order(apply(Rat_tpm_female.revise[,1:5], 1, mean)),]
select = apply(Rat_tpm_female.revise, 1, function(x) sum( x == 0 ) == 0 )
Rat_tpm_female.revise = Rat_tpm_female.revise[select,]
k = heatmap.2(as.matrix(log(Rat_tpm_female.revise, base =2 )), trace = "none", col = Colors, scale = "none",Rowv = FALSE, Colv = FALSE, main = "Female 6Mo vs 6Wo heatmap [Age]")
Rat_tpm_female.revise.norm = t(as.data.frame(k$carpet))
Rat_tpm_female.revise.norm = Rat_tpm_female.revise.norm[order(apply(Rat_tpm_female.revise.norm[,6:10] ,1 ,mean)),]
heatmap.2(as.matrix(Rat_tpm_female.revise.norm), trace = "none", col = Colors, scale = "none", Rowv = FALSE, Colv = FALSE, main = "Female 6Mo vs 6Wo heatmap [log2FC >")
