library(ggplot2)
library(reshape2)
tbl = read.table('SH-SY5Y_IronToxicityTest_FeCl2_24h_WST-1_results.csv', row.names = 1, header = TRUE)
conc.list = c('Control','100uM','200uM','500uM','1mM','2mM','5mM','10mM')
colnames(tbl) = conc.list
tbl = apply(t(tbl),1,function(x){(x/tbl[,1])*100})
val.Mean = apply(tbl,2,mean)
val.Mean = c(as.numeric(val.Mean))
val.sd = apply(tbl,2,sd)
val.sd = c(as.numeric(val.sd))
tbl.dot = tbl
#Melt table for ggplot2 plotting
tbl.dotplot = melt(tbl.dot)
#Synchronize plot name to dotplot variable name
tbl.new = t(data.frame(val.Mean))
tbl.plot = melt(tbl.new)
tbl.plot$variable = conc.list
tbl.plot$sd = val.sd
tbl.plot = tbl.plot[,3:5]
colnames(tbl.plot) = c('variable', 'value', 'sd')
colnames(tbl.dotplot)[2] = c('variable')
tbl.plot$value = factor(tbl.plot$value, levels = tbl.plot$value)
p = ggplot(tbl.plot,aes(x = value, y = variable)) + geom_bar(stat = 'identity', position = 'dodge', width = 0.5, colour = 'black',
                                                             fill = 'white')+
  geom_errorbar(aes(x = value, ymin = variable - sd, ymax = variable + sd), position = position_dodge(0.9), width = 0.15)+
  geom_point(aes(x = variable, y = value),tbl.dotplot, position = position_dodge(0.9))+
  theme(text = element_text(size = 20), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15), axis.text.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank()) + 
  ylab('% WST-1 Absorbance\n(compared to control)') + scale_y_continuous(expand = expand_scale(mult = c(0, .1)))

plot(p)
