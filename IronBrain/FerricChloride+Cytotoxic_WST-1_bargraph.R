library(ggplot2)
library(reshape2)
library(grid)
tbl = read.table('SH-SY5Y_IronToxicityTest_FeCl2_24h_WST-1_results.csv', row.names = 1, header = TRUE)
conc.list = c('Control','100','200','500','1','2','5','10')
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
pdf(file = "Figure4_1_SH-SY5Y_IronTreat_WST-1_Cytotoxicity_barplot_v1.pdf", width = 3.1, height = 3)
p = ggplot(tbl.plot,aes(x = value, y = variable)) + geom_bar(stat = 'identity', position = 'dodge', width = 0.5, colour = 'black',
                                                             fill = 'white', lwd = 0.8)+
  geom_errorbar(aes(x = value, ymin = variable - sd, ymax = variable + sd), position = position_dodge(0.9), width = 0.3,size = 0.8)+
  geom_point(aes(x = variable, y = value),tbl.dotplot, position = position_dodge(0.9), size = 0.5, fill = 'white', shape = 21)+
  theme(text = element_text(size = 9), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 9), axis.text.y = element_text(colour = 'black', size = 9),
        axis.title.x = element_text(colour = 'black', size = 9), axis.ticks = element_line(size = 0.9)) + 
  ylab('% WST-1 Absorbance\n(compared to control)') + xlab('FeCl2 Concentration (mM)')+scale_y_continuous(expand = expand_scale(mult = c(0, .1)))

plot(p)

dev.off()
