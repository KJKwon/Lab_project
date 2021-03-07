library(ggplot2)
tbl = read.table('20210306_GADD45B_siRNA_test_results.csv',  header = TRUE, sep = '\t')
Avg.mean = c()
Avg.sd = c()
Avg.group = character()
Avg.Target = character()
for(i in 0:(length(tbl$Cq)/3-1)){
  temp.mean = mean(tbl$Cq[(i*3+1):((i+1)*3)])
  temp.sd = sd(tbl$Cq[(i*3+1):((i+1)*3)])
  temp.group = as.character(tbl$Group[(i*3+1)])
  temp.taget = as.character(tbl$Target[(i*3+1)])
  Avg.mean = c(Avg.mean,temp.mean)
  Avg.sd = c(Avg.sd,temp.sd)
  Avg.group = c(Avg.group,temp.group)
  Avg.Target = c(Avg.Target,temp.taget)
}
tbl.mean = data.frame(Avg.mean,Avg.sd,Avg.group,Avg.Target)

p = ggplot(tbl.mean,aes(x = Avg.group, y = Avg.mean, fill = Avg.Target))+
  geom_bar(position = position_dodge(),stat = 'identity')+
  geom_point(aes(x = Group, y = Cq, color = Target),
             data = tbl,stat='identity', position = position_dodge(0.9), inherit.aes = FALSE,
             size = 1)+
  geom_errorbar(aes(ymin = Avg.mean - Avg.sd, ymax = Avg.mean + Avg.sd), position =position_dodge(0.9), width = 0.15)+
  scale_color_manual(values = c('black','black'))+
    theme(legend.title = element_blank(), plot.title = element_text(size = 30, hjust = 0.5,margin = margin(0,0,40,0)), 
        panel.background = element_blank(), axis.line = element_line(size = 1, colour = 'black'), 
        legend.position = 'none', axis.text.x = element_text(size = 15, colour = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15, colour = 'black'), axis.title.y = element_text(size = 20))+ 
  xlab(NULL) +ylab('Cq')+ 
  scale_y_continuous(expand = c(0,0), limits = c(0, 25))+
  labs(title = "Gadd45b siRNA 48h (n = 2)")+ scale_fill_grey()
plot(p)
