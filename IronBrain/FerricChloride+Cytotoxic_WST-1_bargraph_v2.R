library(ggplot2)
library(reshape2)
library(RColorBrewer)

#Process multiple table at once
Absorb.tbl.nm.total <- list.files(path = '.', pattern = '_SelectedAbsorbance.txt$')
Absorb.tbl.list <- list()
Absorb.samples <- c()
tmp.dotplot = data.frame()
for(tmp.file.nm in Absorb.tbl.nm.total){
  tmp.tbl <- read.table(tmp.file.nm, sep = '\t', row.names = 1, header = TRUE)
  ##Relative value option
  #tmp.tbl = (tmp.tbl/tmp.tbl[,1])*100
  Absorb.tbl.list <- append(Absorb.tbl.list,list(tmp.tbl))
  Absorb.samples <- c(Absorb.samples,unlist(strsplit(tmp.file.nm,'_'))[2])
  tmp.dotplot = rbind(tmp.dotplot,tmp.tbl)
}

#table of mean value for bar height
val.mean = lapply(Absorb.tbl.list, FUN = function(x) apply(x, 2, mean))
val.sd = lapply(Absorb.tbl.list, FUN = function(x) apply(x, 2, sd))
val.mean = data.frame(matrix(unlist(val.mean), nrow=length(val.mean), byrow=T))
colnames(val.mean) = c('Control','500uM','2mM')
val.mean$group = factor(Absorb.samples, levels = c('NonTreat','NegativeControl','Ptp4a3','Vegfb'))
val.mean = melt(val.mean)

#table of standard deviation for error bar
val.sd = data.frame(matrix(unlist(val.sd), nrow=length(val.sd), byrow=T))
colnames(val.sd) = c('Control','500uM','2mM')
val.sd$group = factor(Absorb.samples, levels = c('NonTreat','NegativeControl','Ptp4a3','Vegfb'))
val.sd = melt(val.sd)

#table of each point
colnames(tmp.dotplot) = c('Control', '500uM', '2mM')
tmp.dotplot = melt(tmp.dotplot)
tmp.dotplot$group = factor(rep(rep(Absorb.samples,each = 5),3), levels = c('NonTreat','NegativeControl','Ptp4a3','Vegfb'))
tmp.dotplot$value = as.numeric(as.character(tmp.dotplot$value))

##Relative value option
#val.mean = val.mean[val.mean$variable != 'Control',]
#val.sd = val.sd[val.sd$variable != 'Control',]
#tmp.dotplot = tmp.dotplot[tmp.dotplot$variable != 'Control',  ]

##Final input table
val.total = val.mean
val.total$sd = val.sd$value

##plotting
p = ggplot(val.total, aes(x = variable, y = value, fill = group)) + geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = brewer.pal(4,"RdYlBu")) + geom_errorbar(aes(ymin = value - sd, ymax = value + sd, width = 0.3),
                                                                     position = position_dodge(width = 0.9)) +
  geom_point(aes(x = variable, y = value, group = group), tmp.dotplot, position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.3))+
  theme(text = element_text(size = 12), legend.title = element_blank(), panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 1, colour = 'black'),
        axis.text.x = element_text(colour = 'black', size = 12), axis.text.y = element_text(colour = 'black', size = 9),
        axis.title.x = element_text(colour = 'black', size = 12), axis.ticks = element_line(size = 0.9),
        plot.title = element_text(hjust = 0.9, size  = 15)) + 
  ggtitle('siRNA treated SH-SY5Y WST-1 Relative absorbance ')+
  ylab('Relative WST-1 absorbance\n(Compared to Control %)') + xlab(expression('FeCl'[2]*' Concentration'))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
plot(p)
