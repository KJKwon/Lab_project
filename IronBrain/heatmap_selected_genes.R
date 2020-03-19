library(gplots)
library(RColorBrewer)
tbl = read.table('B_Up_C_Down_edgeR_QLF_clean.GeneSymbol.txt', header = TRUE, row.names = 1)
tbl.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneID.GeneSymbol.txt', header = TRUE)
tbl.candi = tbl.rpkm[tbl.rpkm$GeneName %in% rownames(tbl),]
rownames(tbl.candi) = tbl.candi$GeneName
tbl.candi = tbl.candi[,c(-1)]
tbl.candi = tbl.candi[rownames(tbl[order(tbl$logFC),]),]
#tbl.candi.sorted = tbl.candi[order(apply(tbl.candi,1,mean)),]
#tbl.candi.sorted = tbl.candi.sorted[apply(log(tbl.candi.sorted) < 8 , 1, all),]
tmp_cor = dist(1-cor(tbl.candi, method = 'spearman'))
tmp_clust = hclust(tmp_cor, method = 'average')
tmp_plot = as.dendrogram(tmp_clust)
#gene_cor = dist(1-cor(t(tbl.candi), method = 'spearman'))
#gene_clust = hclust(gene_cor, method = 'average')
#gene_plot = as.dendrogram(gene_clust)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)
heatmap.2(as.matrix(tbl.candi), Colv = tmp_plot, trace = "none",col = Colors, , scale = 'row', Rowv = FALSE)
par(mar = c(1,1,1,3))
heatmap.2(as.matrix(tbl.candi), trace = "none", col = Colors, scale = 'row', Rowv = FALSE, Colv = FALSE, 
          main = "Female 6Mo vs 6Wo heatmap [Age]", keysize = 1,lhei = c(1,5), lwid = c(1.2,5),margins = c(5.5,7))
