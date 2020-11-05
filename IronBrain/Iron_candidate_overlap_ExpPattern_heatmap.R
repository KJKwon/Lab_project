library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(dendextend)
tbl.gene = read.table('Iron_responsive_candidate_table_6Movs6Wo+IronTreat_MGISeq_clean_One2One_concordant.txt', header = FALSE, stringsAsFactors = FALSE, sep = '\t')
tbl.human.rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneName.txt', header = TRUE, row.names = 1)
colnames(tbl.human.rpkm) = c('10mM_1','10mM_2','10mM_3','2mM_1','2mM_2','2mM_3',
                             '500uM_1','500uM_2','500uM_3','Control_1','Control_2','Control_3')
tbl.human.rpkm = tbl.human.rpkm[,c(10:12,7:9,4:6,1:3)]
selected.row = match(tbl.gene$V2, rownames(tbl.human.rpkm))
tbl.human.rpkm.selected = tbl.human.rpkm[selected.row,]
tbl.human.rpkm.scaled = as.data.frame(t(apply(tbl.human.rpkm.selected, 1, scale)))
colnames(tbl.human.rpkm.scaled) = colnames(tbl.human.rpkm)
breaksList = seq(-3.5,3.5, by = 0.005)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(length(breaksList)))
#pheatmap dendrogram is more clear to analyze correlation
ph <- pheatmap(as.matrix(tbl.human.rpkm.scaled), scale = 'none', cluster_rows = TRUE, cluster_cols = FALSE,
               fontsize = 8, border_color = NA, show_rownames = TRUE,
               width = 3, height = 4,treeheight_row = 2.5)
tmp_clust = as.dendrogram(ph$tree_row)
tmp_clust = set(tmp_clust, "branches_lwd",2)
col_fun <- colorRamp2(seq(-3.5,3.5, by = 0.005), rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(length(breaksList))))
pdf(file = "Figure5_2_RatBrain_6Wovs6Mo_MGISeq_clean_One2One_RatHeatmap.pdf", width = 4.9, height = 2.3)
p = Heatmap(as.matrix(tbl.rat.rpkm.scaled), cluster_rows = tmp_clust, cluster_columns = FALSE, col = col_fun,
            heatmap_legend_param = list(title = "",labels_gp = gpar(fontsize = 8),grid_width = unit(0.3,"cm"),
                                        legend_height = unit(3, "cm")), 
            row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
            row_dend_width = unit(1,'cm'), show_row_names = TRUE, show_row_dend = FALSE,
            height = unit(1.5, 'inch'))
tmp_clust = row_dend(p)
tmp_clust = set(tmp_clust, "branches_lwd",2)
draw(p)
dev.off()

  tbl.rat.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', header = TRUE,row.names = 1)
colnames(tbl.rat.rpkm) = c(gsub('^Rat_15M_','M-15Mo-',colnames(tbl.rat.rpkm)[1:8]),gsub('^Rat_6M_','M-6Mo-',colnames(tbl.rat.rpkm)[9:16]),
                           gsub('^Rat_6M_','F-6Mo-',colnames(tbl.rat.rpkm)[17:26]),gsub('^Rat_6W_','F-6Wo-',colnames(tbl.rat.rpkm)[27:36]))
#tbl.rat.rpkm.selected = tbl.rat.rpkm[match(tbl.gene$V6, rownames(tbl.rat.rpkm)),]
tbl.rat.rpkm.selected = tbl.rat.rpkm[match(tbl.gene$V6, rownames(tbl.rat.rpkm)),]
tbl.rat.rpkm.scaled = as.data.frame(t(apply(tbl.rat.rpkm.selected, 1, scale)))
colnames(tbl.rat.rpkm.scaled) = colnames(tbl.rat.rpkm.selected)
tbl.rat.rpkm.scaled = tbl.rat.rpkm.scaled[,c(36:1)]
#breaksList = seq(-3.5,3.5, by = 0.005)
