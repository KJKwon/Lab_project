library(ggplot2)
library(ggrepel)
tbl = read.table('DAY8vsTPP-GA_DAY8_DE_output.csv', sep = '\t', header = TRUE)
cut_FC = 1.00
cut_FDR = 0.05
down_adipokine = c('Cfd','Retn','Adipoq','Rarres2','Igf1','Tnfrsf21','Igf2','Tnfaip6','Lcn2','Cxcl10','Nampt','Vegfc')
up_adipokine = c('Tnfrsf12a','Ccl5','Ctsl')
plain_adipokine = c('Ddit3')
tbl_ready = tbl[,-c(3,4)]
tbl_ready$group = 'non_adipokine'
tbl_ready$group[tbl_ready$gene_list %in% up_adipokine] = 'up_adipokine'
tbl_ready$group[tbl_ready$gene_list %in% down_adipokine] = 'down_adipokine'
tbl_ready$group[tbl_ready$gene_list %in% plain_adipokine] = 'plain_adipokine'
colour = ifelse(tbl_ready$FDR >= 0.05 | abs(tbl_ready$logFC) < 1, '#CCCCCC', ifelse(tbl_ready$logFC > 1, '#FF3300', '#0066FF'))
tbl_ready = cbind(tbl_ready, colour)
tbl_ready_ordered = tbl_ready[order(tbl_ready$FDR),]
tbl_ready_ordered$gene_list = factor(tbl_ready_ordered$gene_list, levels(tbl_ready_ordered$gene_list)) 
p = ggplot(tbl_ready_ordered, aes(x = logFC, y = -log10(FDR), color = colour,label = gene_list)) + geom_point(size = 1) + scale_color_manual(values = c('#FF0000','#CCCCCC', '#0066FF'))+
  theme_linedraw() + theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlim(-10,5) + geom_vline(xintercept = c(1,-1), linetype = "twodash") + geom_hline(yintercept = -log10(0.05), linetype = "twodash")+
  geom_point(aes(x = logFC, y = -log10(FDR)),subset(tbl_ready,group == 'up_adipokine'), colour = '#0066FF', size = 1)+
  geom_point(aes(x = logFC, y = -log10(FDR)),subset(tbl_ready,group == 'down_adipokine', colour = '#FF3300', size = 1))+
  geom_text_repel(data = subset(tbl_ready, group == 'up_adipokine'), colour = 'black',
                  nudge_x = c(0.1,2,1), nudge_y = c(30,5,-10), fontface = 'bold')+
  geom_text_repel(data = subset(tbl_ready, group == 'down_adipokine'), colour = 'black',
                  nudge_x = c(0.5,0.5,0.5,-1,-1,-2,-3,-3,-0.7,-3,-3,-1.5), nudge_y=c(0,0,0,1.5,1,1,15,8,50,5,-15,-12), fontface = 'bold')+
  geom_text_repel(data = subset(tbl_ready, group == 'plain_adipokine'), colour = 'black', nudge_x = -0.8 , nudge_y = 10)
plot(p)
#  geom_text(aes(label = ifelse(group == 'adipokine' , as.character(gene_list), '')))
#nudge_y = 45+log10(subset(tbl_ready, group == 'up_adipokine')$FDR)
