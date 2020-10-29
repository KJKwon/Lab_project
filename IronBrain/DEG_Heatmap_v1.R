library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dendextend)
#HUMAN
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
breaksList = seq(-3.5,3.5, by = 0.005)
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(length(breaksList)))

#RAT
tbl.gene = read.table('6Month_vs_6Week.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.glmQLFit.Robust.edgeR_geneID_clean.txt', header = TRUE,
                      stringsAsFactors = FALSE, sep = '\t')
tbl.rat.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', header = TRUE,row.names = 1)
colnames(tbl.rat.rpkm) = c(gsub('^Rat_15M_','M-15Mo-',colnames(tbl.rat.rpkm)[1:8]),gsub('^Rat_6M_','M-6Mo-',colnames(tbl.rat.rpkm)[9:16]),
                           gsub('^Rat_6M_','F-6Mo-',colnames(tbl.rat.rpkm)[17:26]),gsub('^Rat_6W_','F-6Wo-',colnames(tbl.rat.rpkm)[27:36]))
#tbl.rat.rpkm.selected = tbl.rat.rpkm[match(tbl.gene$V6, rownames(tbl.rat.rpkm)),]
tbl.rat.rpkm.selected = tbl.rat.rpkm[match(rownames(tbl.gene), rownames(tbl.rat.rpkm)),]
tbl.rat.rpkm.scaled = as.data.frame(t(apply(tbl.rat.rpkm.selected, 1, scale)))
colnames(tbl.rat.rpkm.scaled) = colnames(tbl.rat.rpkm.selected)
tbl.rat.rpkm.scaled = tbl.rat.rpkm.scaled[,c(36:1)]
#breaksList = seq(-3.5,3.5, by = 0.005)
tmp_clust = row_dend(p)
tmp_clust = set(tmp_clust, "branches_lwd",2)
col_fun <- colorRamp2(c(-3,0,3), brewer.pal(n=3, name = "RdBu"))

pdf(file = "Figure3_5_RatBrain_6Wo6Mo15Mo_DEG_heatmap.pdf", width = 4.6, height = 4)
p = Heatmap(as.matrix(tbl.rat.rpkm.scaled), cluster_rows = tmp_clust, cluster_columns = FALSE, col = col_fun, 
            heatmap_legend_param = list(title = "",labels_gp = gpar(fontsize = 8),grid_width = unit(0.3,"cm"),
                                        legend_height = unit(2, "cm")), 
            row_names_gp = gpar(fontsize = 8.2), column_names_gp = gpar(fontsize = 8),
            column_dend_height = unit(0.2,"cm"), row_dend_width = unit(1,'cm'), show_row_names = FALSE)
p
dev.off()
