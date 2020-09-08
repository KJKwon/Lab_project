library(ggplot2)
library(reshape2)
tbl = read.table('SH-SY5Y_IronTreat_rpkm_GeneName.txt', row.names = 1, header = TRUE)
tbl = tbl[,c(2,3,4,5,7,8)]
Gene.list = c('GPX4','SLC3A2','SLC7A11','NFE2L2')
tbl.select = tbl[rownames(tbl) %in% Gene.list, ]
tbl.select = log(tbl.select, base = 2)
tbl.select = cbind(GeneName = rownames(tbl.select), tbl.select)
empty.select = data.frame('Control_NA' = c(rep(NA,length(Gene.list))),
                          'SH.SY5Y.2mM_NA' = c(rep(NA,length(Gene.list))),
                          'SH.SY5Y.10mM_NA' = c(rep(NA,length(Gene.list))))
tbl.select = cbind(tbl.select,empty.select)
tbl.barplot = melt(tbl.select, id.vars = 'GeneName')
tbl.barplot = tbl.barplot[order(match(tbl.barplot$variable,c(names(colour.var)))),]
tbl.barplot$variable = factor(tbl.barplot$variable, levels = c(names(colour.var)))
colour.var = c('Control_2' = '#F8766D','Control_3' = '#ED8141', 'Control_NA' = "#CCCCCC",
               'SH.SY5Y.2mM_1' = '#00B81F','SH.SY5Y.2mM_2' = '#00BC59', 'SH.SY5Y.2mM_NA' = '#CCCCCC',
               'SH.SY5Y.10mM_1' = '#00C0B0','SH.SY5Y.10mM_2' = '#00BDD0', 'SH.SY5Y.10mM_NA' = '#CCCCCC')
p = ggplot(tbl.barplot,aes(x = GeneName, y = value, fill = variable)) + 
  geom_bar(stat='identity', position = 'dodge', color = 'black') +
  scale_fill_manual(values = colour.var, breaks = colnames(tbl))+
  theme(text= element_text(size = 15), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  labs(x = NULL, y = expression(log[2](RPKM)))+
  scale_y_continuous(expand = expand_scale(mult = c(.05, .1)))
plot(p)
