library(RColorBrewer)
library(pheatmap)
tbl.candi = read.table('Iron_responsive_candidate_table_6Movs6Wo_clean+IronTreat_MGISeq_One2One_concordant.txt', header = TRUE,
                       stringsAsFactors = FALSE, sep = '\t')
#tbl.candi = tbl.candi[tbl.candi$V7 >0,]
tbl.human.rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneName.txt', header = TRUE, row.names = 1)
colnames(tbl.human.rpkm) = c('SH-SY5Y_10mM_1','SH-SY5Y_10mM_2','SH-SY5Y_10mM_3','SH-SY5Y_2mM_1','SH-SY5Y_2mM_2','SH-SY5Y_2mM_3',
                             'SH-SY5Y_500uM_1','SH-SY5Y_500uM_2','SH-SY5Y_500uM_3','Control_1','Control_2','Control_3')
tbl.human.rpkm.scaled = as.data.frame(t(apply(tbl.human.rpkm, 1, scale)))
colnames(tbl.human.rpkm.scaled) = colnames(tbl.human.rpkm)
tbl.human.rpkm.scaled = tbl.human.rpkm.scaled[,c(10:12,7:9,4:6,1:3)]
tbl.human.rpkm.selected = tbl.human.rpkm.scaled[match(tbl.candi$GeneName, rownames(tbl.human.rpkm)),]
tbl.rat.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', header = TRUE,row.names = 1)
tbl.rat.rpkm.scaled = as.data.frame(t(apply(tbl.rat.rpkm.selected, 1, scale)))
tbl.rat.rpkm.selected = tbl.rat.rpkm[match(tbl.candi$V6, rownames(tbl.rat.rpkm)),]
colnames(tbl.rat.rpkm.scaled) = colnames(tbl.rat.rpkm.selected)
tbl.rat.rpkm.scaled = tbl.rat.rpkm.scaled[,c(36:1)]
breaksList = seq(-3.5,3.5, by = 0.005)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(length(breaksList)))
p = pheatmap(as.matrix(tbl.human.rpkm.selected), scale = 'none', cluster_rows = TRUE, cluster_cols = TRUE,
         color = col, breaks = breaksList, fontsize = 12, border_color = NA, cellwidth = 15, cellheight = 12, show_rownames = TRUE)
temp_col = p$tree_col
temp_col$order = as.integer(c(1,2,3,4,5,6,7,8,9,12,10,11))
p = pheatmap(as.matrix(tbl.human.rpkm.selected), scale = 'none', cluster_rows = TRUE, cluster_cols = temp_col,
             color = col, breaks = breaksList, fontsize = 12, border_color = NA, cellwidth = 15, cellheight = 12, show_rownames = TRUE)
