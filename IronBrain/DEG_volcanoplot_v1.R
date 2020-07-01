library(ggplot2)
tbl = read.table('SH-SY5Y_IronTreat_2mMvs10mM.3pairs.edgeRQLF_output.txt', sep = '\t', header = TRUE)
tbl.sig = read.table('SH-SY5Y_IronTreat_2mMvs10mM.3pairs.edgeRQLF_output_clean.txt', sep = '\t', header = TRUE)
cut_FC = 1.00
cut_FDR = 0.05
IronTreat.up = rownames(tbl.sig[tbl.sig$logFC >= 1,])
IronTreat.down = rownames(tbl.sig[tbl.sig$logFC <= -1,])
tbl = cbind(gene_list = rownames(tbl), tbl)
rownames(tbl) = c()
tbl_ready = tbl[,-c(3,4,5)]
tbl_ready$group = NA
tbl_ready$group[tbl_ready$gene_list %in% IronTreat.down] = '2mM treat up'
tbl_ready$group[tbl_ready$gene_list %in% IronTreat.up] = '10mM treat up'
tbl_ready_ordered = tbl_ready[order(tbl_ready$group),]
p = ggplot(tbl_ready_ordered, aes(x = logFC, y = -log10(FDR), colour = group))+ geom_point(size = 1) +
  scale_color_manual(breaks = c('2mM treat up','10mM treat up'), values = c('#FF3300', '#0066FF') )+
  geom_point(aes(x = logFC, y = -log10(FDR)),subset(tbl_ready_ordered,is.na(tbl_ready_ordered$group)), colour = '#CCCCCC', size = 1)+
  theme_linedraw() + theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = c(1,-1), linetype = "twodash") + geom_hline(yintercept = -log10(0.05), linetype = "twodash")+
  scale_x_continuous(breaks = c(-7:7),limits = c(-7,7))+ labs(x = expression(logFC), y = expression(-log[10](FDR)))+
  ggtitle("2mM vs 10mM Iron treat") + theme(plot.title = element_text(size= 20, face = 'bold'),
                                                  axis.title.x = element_text(size = 15, face = 'bold'),
                                                  axis.text.x = element_text(size = 13),
                                                  axis.title.y = element_text(size = 15, face = 'bold'),
                                                  axis.text.y = element_text(size = 13))+
  theme(legend.position = c(0.88,0.93), legend.title = element_blank(), legend.text = element_text(size = 12))
plot(p)
