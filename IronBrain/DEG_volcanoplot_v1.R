library(ggplot2)
tbl = read.table('SH-SY5Y_IronTreat_Controlvs10mM_MGISeq.edgeRQLF_Robust_GeneName_output.txt', sep = '\t', header = TRUE)
tbl.sig = read.table('SH-SY5Y_IronTreat_Controlvs10mM.MGISeq.edgeRQLF_Robust_GeneName_output_clean.txt', sep = '\t', header = TRUE)
cut_FC = 1.00
cut_FDR = 0.05
IronTreat.up = rownames(tbl.sig[tbl.sig$logFC >= 1,])
IronTreat.down = rownames(tbl.sig[tbl.sig$logFC <= -1,])
tbl = cbind(gene_list = rownames(tbl), tbl)
rownames(tbl) = c()
tbl_ready = tbl[,-c(3,4,5)]
tbl_ready$group = NA
tbl_ready$group[tbl_ready$gene_list %in% IronTreat.down] = 'Control up'
tbl_ready$group[tbl_ready$gene_list %in% IronTreat.up] = '10mM treat up'
tbl_ready_ordered = tbl_ready[order(tbl_ready$group),]
#tmp = as.data.frame(t(c('ENSG00000000000', -7.00, 0.000199, '10mM treat up')))
#colnames(tmp) = colnames(tbl_ready_ordered)
#tbl_ready_ordered = rbind(tbl_ready_ordered, tmp)
pdf(file = "Figure4_5_SH-SY5Y_Control_vs_10mM_Iron_treat_VolcanoPlot_v1.pdf", width = 3.1, height = 3)
p = ggplot(tbl_ready_ordered, aes(x = logFC, y = -log10(FDR), colour = group))+ geom_point(size = 1) +
  scale_color_manual(breaks = c('Control up','10mM treat up'), values = c('#0066FF','#FF3300'))+
  geom_point(aes(x = logFC, y = -log10(FDR)),subset(tbl_ready_ordered,is.na(tbl_ready_ordered$group)), colour = '#CCCCCC', size = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, size = 1))+
  geom_vline(xintercept = c(1,-1), linetype = "twodash", size = 0.8) + geom_hline(yintercept = -log10(0.05), linetype = "twodash",
                                                                                  size = 0.8)+
  scale_x_continuous(breaks = c(-5:5),limits = c(-5,5))+ labs(x = expression(logFC), y = expression(-log[10](FDR)))+
  ggtitle("Control vs 10mM Iron treat") + theme(plot.title = element_text(size= 11, face = 'bold'),
                                            axis.title.x = element_text(size = 9, face = 'bold'),
                                            axis.text.x = element_text(size = 8),
                                            axis.title.y = element_text(size = 9, face = 'bold'),
                                            axis.text.y = element_text(size = 8))+
  theme(legend.position = "none", axis.ticks = element_line(size = 0.8))

plot(p)

dev.off()
