library(RColorBrewer)
library(gplots)
library(pheatmap)
tbl.candi = read.table('Iron_responsive_candidate_table_15Movs6Wo_clean+IronTreat_MGISeq_One2One_concordant.txt', header = FALSE,
                       stringsAsFactors = FALSE)
tbl.human.rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneName.txt', header = TRUE, row.names = 1)
#colnames(tbl.human.rpkm) = c('Control_1','Control_2','Control_3', 'SH-SY5Y_2mM_1','SH-SY5Y_2mM_2','SH-SY5Y_2mM_3','SH-SY5Y_10mM_1','SH-SY5Y_10mM_2','SH-SY5Y_10mM_3')
tbl.rat.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', header = TRUE,row.names = 1)
tbl.human.rpkm.selected = tbl.human.rpkm[match(tbl.candi$V1, rownames(tbl.human.rpkm)),]
tbl.rat.rpkm.selected = tbl.rat.rpkm[match(tbl.candi$V4, rownames(tbl.rat.rpkm)),c(1:8,27:36)]
tbl.rat.rpkm.selected = tbl.rat.rpkm.selected[apply(tbl.rat.rpkm.selected,1, function(x) sum(x >= 1) >= 10),]
tbl.rat.rpkm.scaled = as.data.frame(t(apply(tbl.rat.rpkm.selected, 1, scale)))
colnames(tbl.rat.rpkm.scaled) = colnames(tbl.rat.rpkm.selected)
#rank = apply(tbl.rpkm.scaled.selected[,1:8], 1, function(x) sum(x < 0) ) + apply(tbl.rpkm.scaled.selected[,9:16], 1, function(x) sum(x > 0))
#rank.filtered = names(sort(rank))
tbl.human.rpkm.scaled = as.data.frame(t(apply(tbl.human.rpkm.selected,1, scale)))
colnames(tbl.human.rpkm.scaled) = colnames(tbl.human.rpkm.selected)
rank.filter1 = apply(tbl.human.rpkm.scaled[,1:6], 1, function(x) sum(x > 0 ) >= 5 | sum(x < 0) >= 5) &
  apply(tbl.human.rpkm.scaled[,7:12], 1, function(x) sum(x > 0) >= 5 | sum(x < 0) >= 5)
rank.filter2 = apply(tbl.human.rpkm.scaled[,4:12], 1, function(x) sum(x > 0 ) >= 8 | sum(x < 0) >= 8) &
  apply(tbl.human.rpkm.scaled[,1:3], 1, function(x) sum(x > 0) >= 2 | sum(x < 0) >= 2) 
rank.filter3 = apply(tbl.human.rpkm.scaled[,1:9], 1, function(x) sum(x > 0 ) >= 8 | sum(x < 0) >= 8) &
  apply(tbl.human.rpkm.scaled[,10:12], 1, function(x) sum(x > 0) >= 2 | sum(x < 0) >= 2)
rank.filtered = rownames(tbl.human.rpkm.selected)[rank.filter1 | rank.filter2 | rank.filter3]
manual.MGI.filter = c("GADD45B","VEGFB","CTF1","PTP4A3","SULT1A1","HES5","MPV17L","SFT2D2","EMILIN2","FNDC3B","CDKN1C",
  "HMGCR","GNB4","HMGCS1")
manual.Next.filter = rev(c('HES5','EMILIN2','CPXM1','CDKN1C','DPYSL4','MFAP2','GALNS','EPB41L4B','TMEM192',
                       'XPNPEP3','PTPN11','TIFA'))
tbl.rat.rpkm.scaled = tbl.rat.rpkm.scaled[,c(9:18,1:8)]
tbl.human.rpkm.scaled = tbl.human.rpkm.scaled[,c(10:12,7:9,4:6,1:3)]
tbl.rpkm.scaled.selected = tbl.rat.rpkm.scaled[match(manual.MGI.filter, rownames(tbl.human.rpkm.scaled)),]
tbl.rpkm.scaled.selected = tbl.human.rpkm.scaled[match(rank.filtered, rownames(tbl.human.rpkm.scaled)),]
tbl.rpkm.scaled.selected = tbl.rpkm.scaled.selected[order(apply(tbl.rpkm.scaled.selected[,c(1:3)],1,sum)),]
breaksList = seq(-3.5,3.5, by = 0.005)
col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(length(breaksList))
pheatmap(as.matrix(tbl.rpkm.scaled.selected), scale = 'none', cluster_rows = FALSE, cluster_cols = FALSE,
         color = col, breaks = breaksList, fontsize = 14, border_color = NA, cellwidth = 25,
         cellheight = 14, show_rownames = TRUE)
