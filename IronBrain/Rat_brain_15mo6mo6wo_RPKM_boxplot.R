library(ggplot2)
library(reshape2)
library(RColorBrewer)
#Iron related gene error bar
tbl = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt',sep='\t', header = TRUE)
selected.gene = c('Herpud1')
tbl = tbl[tbl$geneID %in% selected.gene,]
tbl$IronChallenge_Ctrl = apply(tbl[,20:37],1,mean)
tbl$IronChallenge_1mM = apply(tbl[,10:19],1,mean)
tbl$IronChallenge_2mM = apply(tbl[,2:9],1,mean)
#tbl$IronChallenge_3mM = apply(tbl[,11:13],1,mean)
IronTreat_Ctrl_se = apply(tbl[,20:37],1,function(x){sd(x)/sqrt(10)})
IronTreat_1mM_se = apply(tbl[,10:19],1,function(x){sd(x)/sqrt(18)})
IronTreat_2mM_se = apply(tbl[,2:9],1,function(x){sd(x)/sqrt(8)})
#IronTreat_3mM_se = apply(tbl[,11:13],1,function(x){sd(x)/sqrt(3)})
Iron_se = c(IronTreat_Ctrl_se, IronTreat_1mM_se, IronTreat_2mM_se)
tbl.dot = tbl[,-c(38:40)]
tbl.dotplot = melt(tbl.dot)
tbl.new = tbl[,c(1,38:40)]
tbl.barplot = melt(tbl.new)
tbl.barplot$variable = factor(c('6wo','6mo', '15mo'), 
                              levels = c('6wo','6mo', '15mo'))
tbl.barplot$se = Iron_se
tbl.dotplot$variable = factor(rep(c('15mo','6mo','6wo'),c(8,18,10)),
                         levels = c('6wo','6mo','15mo'))
#yticks = c(1,50,100,150,200,250,500,1000,2000,3000,4000)
#trans= function(x){pmin(x,250)+0.1*pmax(x-250,0)}
p = ggplot(tbl.dotplot,aes(x = variable, y = value)) + geom_boxplot(aes(fill = variable))+
  geom_point( position = position_dodge(0.9))+
  theme(legend.title = element_blank(), plot.title = element_text(size = 30, hjust = 0.5,margin = margin(0,0,40,0)), 
        panel.background = element_blank(), axis.line = element_line(size = 1, colour = 'black'), 
        legend.position = 'none', axis.text.x = element_text(size = 22, colour = 'black'),
        axis.text.y = element_text(size = 18, colour = 'black'), axis.title.y = element_text(size = 22))+ 
  xlab(NULL) +ylab('RPKM')+ scale_y_continuous(expand = c(0,0), limits = c(0, 70)) + 
  labs(title = "HERPUD1")+ scale_fill_grey()
p
