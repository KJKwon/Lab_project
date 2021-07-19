library(ggplot2)
library(dplyr)
tbl = read.table('210718_IronTreat_Relative_FoldChange.txt', header =TRUE, sep = '\t', row.names = 1)
tbl_dot = melt(tbl)
tbl_dot$group = factor(rep(c('Scramble','CLU','FTL'), each = 15), levels = c('Scramble','CLU','FTL'))
tbl_dot$conc = factor(rep(c('NonTreat','1mM','2mM'), each  = 5, times = 3), levels = c('NonTreat','1mM','2mM'))
tbl_dot <- tbl_dot[-c(29),]
tbl_mean <- melt(apply(tbl, 2, mean))
tbl_mean$value[6] = mean(tbl$CLU_2mM[-c(4)])
tbl_sd <- melt(apply(tbl,2,sd))
tbl_sd$value[6] = sd(tbl$CLU_2mM[-c(4)])
tbl_bar = merge(tbl_mean,tbl_sd, by = 0)
colnames(tbl_bar) <- c('Samples','mean','sd')
tbl_bar$Samples = factor(tbl_bar$Samples, levels = tbl_bar$Samples[c(7,8,9,1,2,3,4,5,6)])
tbl_bar$conc = factor(rep(c('NonTreat','1mM','2mM'), times = 3), levels = c('NonTreat','1mM','2mM'))
tbl_bar$group = factor(rep(c('CLU','FTL','Scramble'), each = 3), levels = c('Scramble','CLU','FTL'))
p = ggplot(tbl_bar, aes(x = group, y = mean, fill = conc)) + geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = brewer.pal(3,"RdYlBu")) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.3),
                                                                     position = position_dodge(width = 0.9)) +
  geom_point(aes(x = group, y = value, group = conc), tbl_dot, position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15))+
  theme(text = element_text(size = 12), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 12), axis.text.y = element_text(colour = 'black', size = 9),
        axis.title.x = element_blank(), axis.ticks = element_line(size = 0.9),
        plot.title = element_text(hjust = 0.5, size  = 15)) + 
  ggtitle('siRNA treated SH-SY5Y WST-1 Relative absorbance ')+
  ylab('Relative WST-1 absorbance\n(Compared to NonTreat %)') +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

plot(p)
