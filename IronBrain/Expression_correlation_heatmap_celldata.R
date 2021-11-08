library(ComplexHeatmap)
library(RColorBrewer)
library(dendextend)
tbl.rpkm = read.table('SH-SY5Y_IronChallenge_4_lineages.human_ens98_longest_bwa_mem_rpkm_GeneName.txt', row.names = 1, header = TRUE)
tbl.rpkm <- tbl.rpkm[,c(1:9)]
total.filter = apply(tbl.rpkm, 1, function(x) sum(x > 1) >= 3)
#male.filter = apply(tbl.rpkm[,1:16], 1, function(x) sum(x >= 1) >= 8)
#female.filter = apply(tbl.rpkm[,17:ncol(tbl.rpkm)], 1, function(x) sum(x >= 1) >= 10)
#total.filter = male.filter & female.filter
tbl.rpkm.clean = tbl.rpkm[total.filter,]
#colnames(tbl.rpkm.clean) = c(gsub('^Rat_15M','M_15Mo',colnames(tbl.rpkm.clean)[1:8]),gsub('^Rat_6M','M_6Mo',colnames(tbl.rpkm.clean)[9:16]),
#gsub('^Rat_6M','F_6Mo',colnames(tbl.rpkm.clean)[17:26]),gsub('^Rat_6W','F_6Wo',colnames(tbl.rpkm.clean)[27:36]))
colnames(tbl.rpkm.clean) = c(gsub('^IronChallenge_','',colnames(tbl.rpkm.clean)))
tbl.rpkm.clean = tbl.rpkm.clean
tbl.rpkm.cor = cor(tbl.rpkm.clean, method = 'spearman')
sample.name = colnames(tbl.rpkm.clean)
sample.group = rep(c('Control', '1mM', '2mM'), c(3,3,3))
#tmp_cor <- dist( 1 - as.matrix(tbl.rpkm.cor))
#tmp_clust <- hclust(tmp_cor, method="average")
#tmp_clust = as.dendrogram(tmp_clust)
#tmp_clust = reorder(tmp_clust,c(36:11,rep(1,10)),agglo.FUN = mean)
tmp_clust = set(new_dend, "branches_lwd",2.5)
sample.group = sample.group[order.dendrogram(new_dend)]
bar.group = HeatmapAnnotation(group = sample.group, annotation_name_side = 'left',
                              col = list(group = c( "Control" = "#E6F598",
                                                   "1mM" = "#ABDDA4", "2mM" = "#3288BD")),
                              annotation_legend_param = list(grid_height = unit(0.3, "cm"), title = "", grid_width = unit(0.3, "cm"),
                                                             row_names_gp = gpar(fontsize = 8), direction = "horizontal"),
                              show_annotation_name = FALSE, height = unit(0.15, "cm"), simple_anno_size_adjust = TRUE,
                              show_legend = FALSE)
ht.main = Heatmap(as.matrix(tbl.rpkm.cor), cluster_rows = tmp_clust,cluster_columns = tmp_clust, top_annotation = bar.group, 
                  heatmap_legend_param = list(legend_height = unit(0.8, "cm"),title = "Spearman correlation",
                                              grid_width = unit(0.8,"cm"),title_gp = gpar(fontsize = 8),
                                              legend_width = unit(5,"cm"),labels_gp = gpar(fontsize = 8),direction = "horizontal"), 
                  row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                  column_dend_height = unit(0.7,"cm"), row_dend_width = unit(0.7,'cm'),
                  show_heatmap_legend = FALSE)

svg(filename = "Figure_4_2_IronSample_cor_plot", width = 3.5, height = 3.5)
ht.main <- draw(ht.main) 
                #merge_legend = TRUE, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()
ht.main
row_order(ht.main)
old_dend <- row_dend(ht.main)
new_dend <- reorder(old_dend, c(5,4,3,2,6,7,10,10,10)) #give weight for orders
