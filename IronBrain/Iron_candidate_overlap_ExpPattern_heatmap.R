library(ComplexHeatmap)
library(RColorBrewer)
tbl.gene = read.table('SH-SY5Y_IronTreat_Controlvs10mM.MGISeq.edgeRQLF_Robust_GeneName_output_clean.txt', header = TRUE,
                      stringsAsFactors = FALSE, sep = '\t')
#tbl.gene = tbl.gene[tbl.gene$V7 >0,]
tbl.human.rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneName.txt', header = TRUE, row.names = 1)
colnames(tbl.human.rpkm) = c('10mM_1','10mM_2','10mM_3','2mM_1','2mM_2','2mM_3',
                             '500uM_1','500uM_2','500uM_3','Control_1','Control_2','Control_3')
tbl.human.rpkm = tbl.human.rpkm[,c(10:12,7:9,4:6,1:3)]
selected.row = match(tbl.gene$V2, rownames(tbl.human.rpkm))
tbl.human.rpkm.selected = tbl.human.rpkm[match(rownames(tbl.gene), rownames(tbl.human.rpkm)),]
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
pdf(file = "Figure4_6_SH-SY5Y_Controlvs10mM_IronTreat_DEG_heatmap.pdf", width = 3.1, height = 3.6)
p = Heatmap(as.matrix(tbl.human.rpkm.scaled), cluster_rows = TRUE, cluster_columns = TRUE, col = col, 
            heatmap_legend_param = list(title = "",at = c(-3,-2,-1,0,1,2,3),labels_gp = gpar(fontsize = 9),grid_width = unit(0.4,"cm"),
                                        labels_gp = gpar(font = 8.1)), 
            row_names_gp = gpar(fontsize = 8.2), column_names_gp = gpar(fontsize = 8.5),
            column_dend_height = unit(0.2,"cm"), row_dend_width = unit(0.5,'cm'))
p
dev.off()
