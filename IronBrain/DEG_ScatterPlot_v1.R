library(ggplot2)
tbl.rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneName.txt', sep = '\t', header = TRUE)
tbl.sig = read.table('SH-SY5Y_IronTreat_Ctrlvs2mM_MGISeq.edgeRQLF_Robust_GeneName_output_clean.txt', sep = '\t', header = TRUE)
cut_FC = 1.00
cut_FDR = 0.05
IronTreat.up = rownames(tbl.sig[tbl.sig$logFC >= 1,])
IronTreat.down = rownames(tbl.sig[tbl.sig$logFC <= -1,])
Iron_10mM = apply(tbl.rpkm[,2:4],1,mean)
Iron_2mM = apply(tbl.rpkm[,5:7],1,mean)
Iron_Control = apply(tbl.rpkm[,11:13],1,mean)
tbl_ready = data.frame(gene_list = tbl.rpkm$geneID,Iron_10mM = Iron_10mM,Iron_2mM = Iron_2mM,Iron_Control = Iron_Control)
tbl_ready$group = 'UnSig'
tbl_ready$group[tbl_ready$gene_list %in% IronTreat.down] = 'Control treat up'
tbl_ready$group[tbl_ready$gene_list %in% IronTreat.up] = '2mM treat up'
tbl_ready_ordered = tbl_ready[order(tbl_ready$group),]
#tmp = as.data.frame(t(c('ENSG00000000000', -7.00, 0.000199, '10mM treat up')))
#colnames(tmp) = colnames(tbl_ready_ordered)
#tbl_ready_ordered = rbind(tbl_ready_ordered, tmp)
pdf(file = "Figure4_7_SH-SY5Y_Control_vs_2mM_Iron_treat_ScatterPlot_v1.pdf", width = 3, height = 2.5)
p = ggplot(tbl_ready) + geom_point(aes(x = log2(Iron_Control+1), y = log2(Iron_2mM+1), colour = group, alpha = 0.1,size = group)) +
  geom_point(data = subset(tbl_ready, group == 'Control treat up'), aes(x = log2(Iron_Control+1), y = log2(Iron_2mM+1), colour = group,size= group))+
  geom_point(data = subset(tbl_ready, group == '2mM treat up'), aes(x = log2(Iron_Control+1), y = log2(Iron_2mM+1), colour = group, size =group))+
  scale_color_manual(breaks = c('Control treat up','2mM treat up','UnSig'), values = c('#0066FF','#FF3300','#CCCCCC'))+
  scale_size_manual(breaks = c('Control treat up','2mM treat up','UnSig'), values = c(1,1,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, size = 1.5))+
  ggtitle("Control vs 2mM") + theme(plot.title = element_text(size= 10, face = 'bold'),
                                            axis.title.x = element_text(size = 9, face = 'bold'),
                                            axis.text.x = element_text(size = 8),
                                            axis.title.y = element_text(size = 9, face = 'bold'),
                                            axis.text.y = element_text(size = 8))+
  theme(legend.position = "none", axis.ticks = element_line(size = 1))

plot(p)

dev.off()
