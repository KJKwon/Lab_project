library(ComplexHeatmap)
library(RColorBrewer)
tbl.rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneName.txt', row.names = 1, header = TRUE)
total.filter = apply(tbl.rpkm[,c(1,2,3,10,11,12)], 1, function(x) sum(x >= 1) >= 3)
#male.filter = apply(tbl.rpkm[,1:16], 1, function(x) sum(x >= 1) >= 8)
#female.filter = apply(tbl.rpkm[,17:ncol(tbl.rpkm)], 1, function(x) sum(x >= 1) >= 10)
#total.filter = male.filter & female.filter
tbl.rpkm.clean = tbl.rpkm[total.filter,]
#colnames(tbl.rpkm.clean) = c(gsub('^Rat_15M','M_15Mo',colnames(tbl.rpkm.clean)[1:8]),gsub('^Rat_6M','M_6Mo',colnames(tbl.rpkm.clean)[9:16]),
#                             gsub('^Rat_6M','F_6Mo',colnames(tbl.rpkm.clean)[17:26]),gsub('^Rat_6W','F_6Wo',colnames(tbl.rpkm.clean)[27:36]))
colnames(tbl.rpkm.clean) = c(gsub('^Fe_Treat_24h_','',colnames(tbl.rpkm.clean)))
tbl.rpkm.clean = tbl.rpkm.clean[,c(10,11,12,7,8,9,4,5,6,1,2,3)]
tbl.rpkm.cor = cor(tbl.rpkm.clean, method = 'spearman')
sample.name = colnames(tbl.rpkm.clean)
sample.group = rep(rev(c('10mM','2mM','500uM','Control')), c(3,3,3,3))
#tbl.sample = data.frame(sample.name,sample.group)
#sample.group = sample.group[row_order(ht.main)]
tmp_cor <- dist( 1 - cor(as.matrix(tbl.rpkm.clean), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
tmp_clust = reorder(tmp_clust, 1:10)
bar.group = HeatmapAnnotation(group = sample.group, annotation_name_side = 'left',
                              col = list(group = c("Control" = "#FFFFCC", "500uM" = "#FED976",
                                                   "2mM" = "#FC4E2A", "10mM" = "#BD0026")),
                              annotation_legend_param = list(grid_height = unit(0.3, "cm"), title = "", grid_width = unit(0.3, "cm"),
                                                             row_names_gp = gpar(fontsize = 9), direction = "horizontal"),
                              show_annotation_name = FALSE, height = unit(0.15, "cm"), simple_anno_size_adjust = TRUE)

ht.main = Heatmap(as.matrix(tbl.rpkm.cor), cluster_rows = tmp_clust,cluster_columns = tmp_clust, top_annotation = bar.group, 
                  heatmap_legend_param = list(legend_height = unit(1.5, "cm"),title = "Spearman correlation",
                                              grid_width = unit(0.2,"cm"),title_gp = gpar(fontsize = 9),
                                              labels_gp = gpar(fontsize = 9),direction = "horizontal"), row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
                  column_dend_height = unit(0.5,"cm"), row_dend_width = unit(0.5,'cm'))
pdf(file = "Figure4_2_SH-SY5Y_IronTreat_correlation_heatmap_v1.pdf", width = 3.1, height = 3.5)
draw(ht.main, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
