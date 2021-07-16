library(ggplot2)
library(reshape2)
library(RColorBrewer)
#Iron related gene error bar
tbl = read.table('SH-SY5Y_IronChallenge_4_lineages.human_ens98_longest_bwa_mem_rpkm_GeneName.txt',sep='\t', header = TRUE)
selected.gene = c('HERPUD1')
tbl = tbl[tbl$geneID %in% selected.gene,]
tbl$IronTreat_ctrl = apply(tbl[,2:4],1,mean)
tbl$IronTreat_1mM = apply(tbl[,5:7],1,mean)
tbl$IronTreat_2mM = apply(tbl[,8:10],1,mean)
tbl$IronTreat_3mM = apply(tbl[,11:13],1,mean)
IronTreat_Ctrl_se = apply(tbl[,2:4],1,function(x){sd(x)/sqrt(3)})
IronTreat_1mM_se = apply(tbl[,5:7],1,function(x){sd(x)/sqrt(3)})
IronTreat_2mM_se = apply(tbl[,8:10],1,function(x){sd(x)/sqrt(3)})
IronTreat_3mM_se = apply(tbl[,11:13],1,function(x){sd(x)/sqrt(3)})
Rat_se = c(IronTreat_Ctrl_se, IronTreat_1mM_se, IronTreat_2mM_se, IronTreat_3mM_se)
tbl.dot = tbl[,-c(14:17)]
tbl.dotplot = melt(tbl.dot)
tbl.new = tbl[,c(1,14:17)]
tbl.barplot = melt(tbl.new)
tbl.barplot$variable = factor(c('Control','1mM','2mM','3mM'), 
                              levels = c('Control','1mM','2mM','3mM'))
tbl.barplot$se = Rat_se
tbl.dotplot$variable = c(rep(c('Control','1mM','2mM','3mM'), each = 3))
#yticks = c(1,50,100,150,200,250,500,1000,2000,3000,4000)
#trans= function(x){pmin(x,250)+0.1*pmax(x-250,0)}
p = ggplot(tbl.barplot,aes(x = variable, y = value)) + geom_bar(aes(fill = variable),
                                                                stat='identity', position = 'dodge', width = 0.5)+
  geom_errorbar(aes(ymin = value - se, ymax = value + se), position =position_dodge(0.9), width = 0.15)+
  geom_point(aes(x = variable, y = value), tbl.dotplot, position = position_dodge(0.9))+
  theme(legend.title = element_blank(), plot.title = element_text(size = 30, hjust = 0.5,margin = margin(0,0,40,0)), 
        panel.background = element_blank(), axis.line = element_line(size = 1, colour = 'black'), 
        legend.position = 'none', axis.text.x = element_text(size = 18, colour = 'black'),
        axis.text.y = element_text(size = 15, colour = 'black'), axis.title.y = element_text(size = 20))+ 
  xlab(NULL) +ylab('RPKM')+ scale_y_continuous(expand = c(0,0), limits = c(0, 100)) + 
  labs(title = "HERPUD1")+ scale_fill_grey()
p
