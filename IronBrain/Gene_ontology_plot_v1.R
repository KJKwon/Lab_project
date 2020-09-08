library(ggplot2)
tbl = read.table('IronTreat_MGISeq_Clean_selected_GO.txt',sep = '\t',header = TRUE)
tbl.clean = tbl[,c(1,3,8)]
colnames(tbl.clean) = c('Pathway', 'Input', 'FDR')
tbl.clean = tbl.clean[order(tbl.clean$FDR, decreasing = TRUE),]
tbl.clean$Pathway = factor(tbl.clean$Pathway, levels = tbl.clean$Pathway)
p = ggplot(tbl.clean,aes(x = Pathway, y = -log(FDR, base = 10), width = 0.5)) +geom_bar(stat = 'identity') + coord_flip() +
  theme_classic()+ scale_y_continuous(limits = c(0,8),expand = c(0,0)) + scale_fill_manual(values = c("#FFCC00"))+
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(), plot.title = element_text(size = 15, face = "bold", hjust = -1))+
  geom_text(aes(label = Input), position = position_dodge(width = 0.9), hjust = -0.5)+
  ylab('-log10(FDR)') + ggtitle("Control vs 10mM DEG (n = 29) selected ontology")
plot(p)
