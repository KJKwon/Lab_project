library(ggplot2)
library(reshape2)
library(minpack.lm)
library(drc)
#Iron related gene error bar
tbl = read.table('SH-SY5Y_IronToxicityTest_FeCl2_24h_WST-1_results.csv',sep='\t',row.names = 1, header = TRUE)
tbl = apply(t(tbl),1,function(x){(x/tbl[,1])*100})
val.Mean = apply(tbl,2,mean)
val.Mean = c(as.numeric(val.Mean))
#Calculate standard error
val.se = apply(tbl,2,function(x){sd(x)/sqrt(length(rownames(tbl)))})
val.se = c(as.numeric(val.se))
tbl.dot = tbl
#Melt table for ggplot2 plotting
tbl.dotplot = melt(tbl.dot)
tbl.dotplot$variable = rep(c(0,100,200,500,1000,2000,5000,10000), times = 1, each = length(rownames(tbl))) 
#Synchronize plot name to dotplot variable name
tbl.new = t(data.frame(val.Mean))
tbl.plot = melt(tbl.new)
tbl.plot$variable = c(0,100,200,500,1000,2000,5000,10000)
tbl.plot$se = val.se
tbl.plot = tbl.plot[,3:5]
colnames(tbl.plot) = c('variable', 'value', 'se')
p = ggplot(tbl.plot,aes(x = log10(value+1), y = variable)) + geom_line() + geom_point(aes(x = log10(variable+1) , y = value), tbl.dotplot)+
  geom_errorbar(aes(ymin = variable - se, ymax = variable + se))+
  theme(text = element_text(size = 15), legend.title = element_blank())+ xlab('log(FeCl2) uM') + ylab('Absorbance(450nm-690nm)')
plot(p)
##Test_set
p = ggplot() + geom_point(aes(x = log10(variable+1), y = value),tbl.dotplot)+ 
#  geom_smooth(data= tbl.dotplot, aes(x = log10(variable+1), y = value), method = 'nls' , 
#              formula = y ~ f(x,xmid,scal), method.args = list(start=c(xmid = 3.4910, scal = -0.6797)),se = FALSE) + 
  geom_smooth(data = tbl.dotplot, aes(x = log10(variable+1), y = value), method = 'auto', se = F)+
  geom_errorbar(aes(x = log10(value+1), ymin = variable - se, ymax = variable + se), tbl.plot)+
  #panel design part
  theme(text = element_text(size = 15), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15), axis.text.y = element_text(colour = 'black', size = 15))+ 
  xlab('log(FeCl2) uM') + ylab('% WST-1 Absorbance\n(compared to control)') +   scale_y_continuous(breaks=c(0,25,50,75,100),
                                                                                                   limits = c(0,100))+
  annotation_logticks(sides = 'b')
plot(p)

##geom_point = scatter_plot, geom_errorbar = error_bar
p = ggplot() + geom_point(aes(x = log10(variable+1), y = value),tbl.dotplot)+ 
  geom_smooth(data= tbl.dotplot, aes(x = log10(variable+1), y = value), method = 'nls' , 
              formula = y ~ SSlogis(x, Asym, xmid, scal), se = FALSE) + 
  geom_errorbar(aes(x = log10(value+1), ymin = variable - se, ymax = variable + se), tbl.plot)+
  #panel design part
  theme(text = element_text(size = 15), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 15), axis.text.y = element_text(colour = 'black', size = 15))+ 
  xlab('log(FeCl2) uM') + ylab('% WST-1 Absorbance\n(compared to control)') +   scale_y_continuous(breaks=c(0,25,50,75,100,125),
                                                                                                   limits = c(0,125))+
  annotation_logticks(sides = 'b')
#scale_y_continuous(breaks = seq(0,125,25))
plot(p)
ggsave(file = "200401_Iron_treat_suction_SH-SY5Y_WST-1_original.pdf", plot = last_plot(), width = 8, height = 6, units = c("in"),
       device = pdf())
dev.off()
