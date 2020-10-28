library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
tbl.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', row.names = 1, header = TRUE)
total.filter = apply(tbl.rpkm[,c(1,2,3,10,11,12)], 1, function(x) sum(x >= 1) >= 3)
#male.filter = apply(tbl.rpkm[,1:16], 1, function(x) sum(x >= 1) >= 8)
#female.filter = apply(tbl.rpkm[,17:ncol(tbl.rpkm)], 1, function(x) sum(x >= 1) >= 10)
#total.filter = male.filter & female.filter
tbl.rpkm.clean = tbl.rpkm[total.filter,]
#colnames(tbl.rpkm.clean) = c(gsub('^Rat_15M','M_15Mo',colnames(tbl.rpkm.clean)[1:8]),gsub('^Rat_6M','M_6Mo',colnames(tbl.rpkm.clean)[9:16]),
                             gsub('^Rat_6M','F_6Mo',colnames(tbl.rpkm.clean)[17:26]),gsub('^Rat_6W','F_6Wo',colnames(tbl.rpkm.clean)[27:36]))
colnames(tbl.rpkm.clean) = c(gsub('^Fe_Treat_24h_','',colnames(tbl.rpkm.clean)))
tbl.rpkm.clean = tbl.rpkm.clean[,c(10,11,12,7,8,9,4,5,6,1,2,3)]
tbl.rpkm.cor = cor(tbl.rpkm.clean, method = 'spearman')
sample.name = colnames(tbl.rpkm.clean)
#sample.group = rep(c('Male 15Mo', 'Male 6Mo', 'Female 6Mo', 'Female 6Wo'), c(8,8,10,10))
tmp_cor <- dist( 1 - as.matrix(tbl.rpkm.cor))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
#tmp_clust = reorder(tmp_clust,c(36:11,rep(1,10)),agglo.FUN = mean)
tmp_clust = set(tmp_clust, "branches_lwd",2.5)
sample.group = sample.group[order.dendrogram(sample.group)]
bar.group = HeatmapAnnotation(group = sample.group, annotation_name_side = 'left',
                              col = list(group = c("Male 15Mo" = "#F46D43", "Male 6Mo" = "#E6F598",
                                                   "Female 6Mo" = "#ABDDA4", "Female 6Wo" = "#3288BD")),
                              annotation_legend_param = list(grid_height = unit(0.3, "cm"), title = "", grid_width = unit(0.3, "cm"),
                                                             row_names_gp = gpar(fontsize = 8), direction = "horizontal"),
                              show_annotation_name = FALSE, height = unit(0.15, "cm"), simple_anno_size_adjust = TRUE)
ht.main = Heatmap(as.matrix(tbl.rpkm.cor), cluster_rows = tmp_clust,cluster_columns = tmp_clust, top_annotation = bar.group, 
                  heatmap_legend_param = list(legend_height = unit(0.8, "cm"),title = "Spearman correlation",
                                              grid_width = unit(0.8,"cm"),title_gp = gpar(fontsize = 8),
                                              legend_width = unit(5,"cm"),labels_gp = gpar(fontsize = 8),direction = "horizontal"), 
                                              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                  column_dend_height = unit(0.7,"cm"), row_dend_width = unit(0.7,'cm'))
pdf(file = "Figure3_2_RatBrain_correlation_heatmap_v1.pdf", width = 5, height = 5.5)
draw(ht.main, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
