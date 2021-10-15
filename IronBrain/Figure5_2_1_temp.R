library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)
library(ggpubr)
tbl.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', row.names = 1, header = TRUE)
candi.gene <- c('Hspa5','Herpud1','Clu')
tbl.rpkm <- tbl.rpkm[,1:26]
tbl.rpkm.clean = tbl.rpkm[match(candi.gene,rownames(tbl.rpkm)),]
tbl.rpkm.clean$Gene <- rownames(tbl.rpkm.clean)
tbl.violin <- melt(tbl.rpkm.clean)
tbl.violin$Age <- rep(c('15M','6M'), times = c(24, 54))
tbl.violin.tmp <- tbl.violin[tbl.violin$Gene == 'Clu',]
p.clu = ggplot(tbl.violin.tmp, aes(x = Age, y = value, fill = Age)) + 
  geom_violin( position = 'dodge',trim = TRUE, scale = 'count') +
  scale_fill_manual(values = brewer.pal(8,"RdYlBu")[c(1,8)])+
  #geom_point(position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.2), show.legend = FALSE,
  #           size = 1)+
  geom_quasirandom(width = 0.1, dodge.width = 0.9, size = 0.5)+
  #geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 0.3, fill = 'black')+
  theme(legend.position = 'none',
        panel.background = element_blank(), panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 9), axis.text.y = element_text(colour = 'black', size = 9),
        axis.title.x = element_blank(), axis.ticks = element_line(size = 0.9),
        axis.title.y = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0,1500))


plot(p.clu)
ggarrange(p.hspa,p.herpud,p.clu, ncol = 3, nrow = 1)
