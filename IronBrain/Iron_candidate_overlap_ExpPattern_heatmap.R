library(RColorBrewer)
library(pheatmap)
tbl.gene = read.table('Iron_responsive_candidate_table_6Movs6Wo_clean+IronTreat_MGISeq_One2One_concordant.txt', header = FALSE,
                       stringsAsFactors = FALSE, sep = '\t')
#tbl.gene = tbl.gene[tbl.gene$V7 >0,]
tbl.human.rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneName.txt', header = TRUE, row.names = 1)
colnames(tbl.human.rpkm) = c('SH-SY5Y_10mM_1','SH-SY5Y_10mM_2','SH-SY5Y_10mM_3','SH-SY5Y_2mM_1','SH-SY5Y_2mM_2','SH-SY5Y_2mM_3',
                             'SH-SY5Y_500uM_1','SH-SY5Y_500uM_2','SH-SY5Y_500uM_3','Control_1','Control_2','Control_3')
tbl.human.rpkm = tbl.human.rpkm[,c(10:12,7:9,4:6,1:3)]
selected.row = match(tbl.gene$V2, rownames(tbl.human.rpkm))
tbl.human.rpkm.selected = tbl.human.rpkm[match(tbl.gene$V2, rownames(tbl.human.rpkm)),]
tbl.human.rpkm.scaled = as.data.frame(t(apply(tbl.human.rpkm.selected, 1, scale)))
colnames(tbl.human.rpkm.scaled) = colnames(tbl.human.rpkm)
tbl.rat.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', header = TRUE,row.names = 1)
colnames(tbl.rat.rpkm) = c(gsub('^Rat_15M','M_15Mo',colnames(tbl.rat.rpkm)[1:8]),gsub('^Rat_6M','M_6Mo',colnames(tbl.rat.rpkm)[9:16]),
                             gsub('^Rat_6M','F_6Mo',colnames(tbl.rat.rpkm)[17:26]),gsub('^Rat_6W','F_6Wo',colnames(tbl.rat.rpkm)[27:36]))
tbl.rat.rpkm.selected = tbl.rat.rpkm[match(tbl.gene$V6, rownames(tbl.rat.rpkm)),]
tbl.rat.rpkm.scaled = as.data.frame(t(apply(tbl.rat.rpkm.selected, 1, scale)))
colnames(tbl.rat.rpkm.scaled) = colnames(tbl.rat.rpkm.selected)
tbl.rat.rpkm.scaled = tbl.rat.rpkm.scaled[,c(36:1)]
breaksList = seq(-3.5,3.5, by = 0.005)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(length(breaksList)))
p = pheatmap(as.matrix(tbl.human.rpkm.scaled), scale = 'none', cluster_rows = TRUE, cluster_cols = FALSE,
             color = col, breaks = breaksList, fontsize = 11, border_color = NA, cellwidth = 15, cellheight = 10, show_rownames = TRUE)
human.order = p$tree_row[["order"]]
tbl.rat.rpkm.scaled.ordered = tbl.rat.rpkm.scaled[human.order,]
temp_col = p$tree_col
temp_col$order = as.integer(c(1,2,3,4,5,6,7,8,9,12,10,11))
p = pheatmap(as.matrix(tbl.rat.rpkm.scaled.ordered), scale = 'none', cluster_rows = FALSE, cluster_cols = FALSE,
             color = col, breaks = breaksList, fontsize = 12, border_color = NA, cellwidth = 15, cellheight = 12, show_rownames = TRUE)
