library(ggplot2)
library(reshape2)
library(minpack.lm)
library(drc)
tbl_Original = read.table('200414_Iron_treat_half_exchange.txt', sep = '\t', row.names = 1, header = TRUE)
tbl_3H_54 = read.table('210206_Iron_treat_half_exchange_3H_54.csv', sep = '\t', row.names = 1, header = TRUE)
tbl_3H_105 = read.table('210206_Iron_treat_half_exchange_3H_105.csv', sep = '\t', row.names = 1, header = TRUE)
tbl_4H_54 = read.table('210206_Iron_treat_half_exchange_4H_54.csv', sep = '\t', row.names = 1, header = TRUE)
tbl_4H_105 = read.table('210206_Iron_treat_half_exchange_4H_105.csv', sep = '\t', row.names = 1, header = TRUE)
list_IronTreat = list(tbl_Original,tbl_3H_54,tbl_3H_105,tbl_4H_54,tbl_4H_105)
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
group =  rep(c('tbl_Original','tbl_3H_54','tbl_3H_105','tbl_4H_54','tbl_4H_105'),8)
val.mean$group = group
temp.lineplot = val.mean
temp.lineplot$group = factor(temp.lineplot$group, levels = c('tbl_Original','tbl_3H_54','tbl_3H_105','tbl_4H_54','tbl_4H_105'))
#temp = val.mean[val.mean$group == 'NT',]
#temp.dotplot = list_IronTreat[[1]]
temp.dotplot <- rbind(tbl_Original,tbl_3H_54,tbl_3H_105,tbl_4H_54,tbl_4H_105)
colnames(temp.dotplot) = c(0,100,200,500,1000,2000,5000,10000)
temp.dotplot = melt(temp.dotplot)
temp.dotplot$group = rep(rep(c('tbl_Original','tbl_3H_54','tbl_3H_105','tbl_4H_54','tbl_4H_105'),each = 5),8)
temp.dotplot$group = factor(temp.dotplot$group, levels = c('tbl_Original','tbl_3H_54','tbl_3H_105','tbl_4H_54','tbl_4H_105'))
temp.dotplot$variable = as.numeric(as.character(temp.dotplot$variable))
p = ggplot(data = temp.lineplot, aes(x = log10(variable+1), y = value, color = group))+
  geom_line(aes(color =group)) + geom_point(aes(x = log10(variable+1), y = value, color = group), temp.dotplot)+
#  geom_errorbar(aes(ymin = value - sd, ymax = value + sd, width = 0.3))+
  scale_color_manual(values=c('Black','Red','Pink','Blue','cyan'), 
                     labels=c('Standard', expression(paste("3H 5X","10"^3)),expression(paste("3H 1X","10"^4)),
                              expression(paste("4H 5X","10"^3)),expression(paste("4H 1X","10"^4))))+
  theme(text = element_text(size = 9), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 9), axis.text.y = element_text(colour = 'black', size = 9),
        axis.title.x = element_text(colour = 'black', size = 9), axis.ticks = element_line(size = 0.9),
        plot.title = element_text(hjust = 0.5, size  = 15)) + 
  ggtitle('WST-1 3H and 4H incubation ')+
  ylab('% WST-1 Absorbance\n(compared to control)') + xlab('FeCl2 Concentration (mM)')+scale_y_continuous(limits = c(-5,100))
plot(p)
