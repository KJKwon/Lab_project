library(ggplot2)
library(scales)
tbl = read.table('6Mo_vs_6Wo_Panther_GeneEnrichmentAnalysis_6Wo_UP_clean.txt',sep = '\t',header = TRUE)
tbl.clean = tbl[,c(1,3,6)]
colnames(tbl.clean) = c('Pathway', 'Input', 'FDR')
tbl.clean = tbl.clean[order(tbl.clean$FDR, decreasing = TRUE),]
tbl.clean$Pathway = c(gsub(pattern = "\\s[(].+[)]",replacement = "",tbl.clean$Pathway))
tbl.clean$Pathway = factor(tbl.clean$Pathway, levels = tbl.clean$Pathway)
pdf(file = "Figure3_6_RatBrain_6Mo_vs_6Wo_GeneEnrichment_plot_v2.pdf", width = 3.3, height = 2)
p = ggplot(tbl.clean,aes(x = Pathway, y = -log(FDR, base = 10), width = 0.5)) +geom_bar(stat = 'identity',width = 1, color = 'black', lwd = 0.8) + 
  coord_flip()+
  theme_classic()+ scale_y_continuous(limits = c(0,8),expand = c(0,0)) + scale_fill_manual(values = c("#FFCC00"))+
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(), plot.title = element_blank(), axis.line = element_line(size = 1, colour = 'black'), 
        axis.ticks.x = element_line(size = 1, colour = 'black'), axis.ticks.y = element_line(size = 1, colour = 'black'))+
  geom_text(aes(label = Input), position = position_dodge(width = 0.5), hjust = -0.5, size = 3)+
  ylab('-log10(FDR)') + ggtitle("6Mo vs 6Wo DEG (n = 29)") +scale_x_discrete(labels = wrap_format(20))

plot(p)
dev.off()
