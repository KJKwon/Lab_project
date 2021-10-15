library(ggplot2)
library(reshape2)
library(RColorBrewer)
tbl_ctrl = read.table('211011_IronTreat_Relative_FoldChange_TimeDependent.csv', header =TRUE, sep = '\t')
tbl_ctrl = tbl_ctrl[tbl_ctrl$Time == '4H',-c(10)]
tbl_treat = read.table('211007_IronTreat_Relative_FoldChange.csv', header = TRUE, sep = '\t')
tbl_treat = tbl_treat[,-c(1)]
tbl_total = rbind(tbl_ctrl, tbl_treat)
tbl_total <- tbl_total[,-c(1,4,7)]
tbl_total$sample <- factor(rep(c('Control','Challenged'),time = c(5,5)), levels = c('Control','Challenged'))
tbl_dot = melt(tbl_total)
tbl_dot$group = factor(rep(c('Scramble','CLU','HERPUD1'), each = 20), levels = c('Scramble','CLU','HERPUD1'))
tbl_dot$conc = factor(rep(c('2mM','5mM'), each  = 10, times = 3), levels = c('2mM','5mM'))
#colnames(tbl_dot) <- c('variable','value','Samples','conc')
#tbl_mean_total = melt(apply(tbl[1:5,1:9], 2, mean))
tbl_mean_ctrl <- melt(apply(tbl_total[1:5,1:6], 2, mean))
tbl_mean_treat <- melt(apply(tbl_total[6:10,1:6], 2, mean))
#tbl_mean_4H <- melt(apply(tbl[11:15,1:6], 2, mean))
tbl_mean_total <- as.data.frame(rbind(tbl_mean_ctrl,tbl_mean_treat))
#tbl_mean_total$group <- factor(rep(c('2H','3H','4H'),each = 6), levels = c('2H','3H','4H'))
#tbl_sd_total = melt(apply(tbl[1:5,1:9], 2, sd))
tbl_sd_ctrl <- melt(apply(tbl_total[1:5,1:6], 2, sd))
tbl_sd_treat <- melt(apply(tbl_total[6:10,1:6], 2, sd))
#tbl_sd_4H <- melt(apply(tbl[11:15,1:6], 2, sd))
tbl_sd_total <- as.data.frame(rbind(tbl_sd_ctrl,tbl_sd_treat))
#tbl_sd_total$group <- factor(rep(c('2H','3H','4H'),each = 6), levels = c('2H','3H','4H'))
tbl_bar = cbind(tbl_mean_total,tbl_sd_total)
#tbl_bar <- tbl_bar[,c(1,3,4)]
tbl_bar$Samples = factor(rep(c('Scramble','CLU','HERPUD1'),each = 2,times = 2 ), levels = c('Scramble','CLU','HERPUD1'))
tbl_bar$conc = factor(rep(c('2mM','5mM'), times = 6), levels = c('2mM','5mM'))
tbl_bar$group = factor(rep(c('Control','Challenged'),time = c(6,6)), levels = c('Control','Challenged'))
colnames(tbl_bar) <- c('mean','sd','group','conc','sample')
p = ggplot(tbl_bar, aes(x = conc, y = mean, fill = group, color = group)) + geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = c('white','grey60','grey40')) + scale_color_manual(values = c('black','black','black'))+
#  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 0.3), position = position_dodge(width = 0.9)) +
  geom_point(aes(x = conc, y = value, group = group), tbl_dot, position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.45),
             show.legend = FALSE, size = 1, shape = 21)+
  theme(text = element_text(size = 15), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 12), axis.text.y = element_text(colour = 'black', size = 9),
        axis.title.x = element_blank(), axis.ticks = element_line(size = 0.6),
        plot.title = element_text(hjust = 0.5, size  = 20, face = 'bold')) + 
  facet_wrap(~sample)+
  ylab('Relative viability\n(compared to NonTreat (%))') +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

plot(p)

test = tbl_dot[tbl_dot$conc == '5mM' & tbl_dot$sample == 'Challenged',]
aov_res <- aov(value~group, data = test)
summary(aov(value~group, data = test))
TukeyHSD(aov_res)
