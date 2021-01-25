library(ggplot2)
library(reshape2)
library(minpack.lm)
library(drc)
tbl_NT = read.table('210119_IronTreat_siRNA_NonTreat_1X104.txt', sep = '\t', row.names = 1, header = TRUE)
tbl_Neg = read.table('210119_IronTreat_siRNA_NegativeControl_1X104.txt', sep = '\t', row.names = 1, header = TRUE)
tbl_Vegfb = read.table('210119_IronTreat_siRNA_Vegfb_1X104.txt', sep = '\t', row.names = 1, header = TRUE)
tbl_Ptp4a3 = read.table('210119_IronTreat_siRNA_Ptp4a3_1X104.txt', sep = '\t', row.names = 1, header = TRUE)
list_IronTreat = list(tbl_NT,tbl_Neg,tbl_Vegfb,tbl_Ptp4a3)
list_IronTreat = lapply(list_IronTreat, FUN = function(x) apply(t(x),1,function(i){(i/x[,1])*100}))
val.mean = lapply(list_IronTreat, FUN = function(x) apply(x, 2, mean))
val.sd = lapply(list_IronTreat, FUN = function(x) apply(x, 2, sd))
val.mean = data.frame(matrix(unlist(val.mean), nrow=length(val.mean), byrow=T))
colnames(val.mean) = c(0,100,200,500,1000,2000,5000,10000)
val.sd = data.frame(matrix(unlist(val.sd), nrow=length(val.sd), byrow=T))
colnames(val.sd) = c(0,100,200,500,1000,2000,5000,10000)
val.mean = melt(val.mean)
val.sd = melt(val.sd)
val.mean$sd = val.sd$value
val.mean$variable = as.numeric(as.character(val.mean$variable))
group =  rep(c('NT','NC','Vegfb','Ptp4a3'),4)
val.mean$group = group
temp = val.mean[val.mean$group == 'NT',]
temp.dotplot = list_IronTreat[[1]]
colnames(temp.dotplot) = c(0,100,200,500,1000,2000,5000,10000)
temp.dotplot = melt(temp.dotplot)
temp.dotplot$Var2 = as.numeric(as.character(temp.dotplot$Var2))
temp.dotplot$group = 'NT'
p = ggplot(data = temp, aes(x = log10(variable+1), y = value, color = group))+
  geom_line(aes(color =group)) + geom_point(aes(x = log10(Var2+1), y = value, color = group), temp.dotplot)+
#  geom_errorbar(aes(ymin = value - sd, ymax = value + sd, width = 0.3))+
  scale_color_manual(values=c('blue'))+
  theme(text = element_text(size = 9), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 9), axis.text.y = element_text(colour = 'black', size = 9),
        axis.title.x = element_text(colour = 'black', size = 9), axis.ticks = element_line(size = 0.9),
        legend.position = 'none', plot.title = element_text(hjust = 0.5, size  = 15)) + 
  ggtitle('NonTreat')+
  ylab('% WST-1 Absorbance\n(compared to control)') + xlab('FeCl2 Concentration (mM)')+scale_y_continuous(limits = c(-5,125))
plot(p)
