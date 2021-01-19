library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(PMCMRplus)
library(lawstat)
library(agricolae)
library(DescTools)
library(ggsignif)
tbl = read.table('20201209_PTP4A3_48h_siRNA_validation_results.txt',  header = TRUE, sep = '\t')
colnames(tbl) = c('Samples','FoldChange')
tbl = tbl[3:6,]
Avg.Negative = mean(tbl$FoldChange[1:2])
Avg.siRNA_48h = mean(tbl$FoldChange[3:4])
#Avg.siRNA.6x = mean(tbl$FoldChange[7:9])
#Avg.siRNA.3x = mean(tbl$FoldChange[10:12])
SD.Negative = sd(tbl$FoldChange[1:2])/sqrt(2)
SD.siRNA_48h = sd(tbl$FoldChange[3:4])/sqrt(2)
#SD.siRNA.6x = sd(tbl$FoldChange[7:9])/sqrt(3)
#SD.siRNA.3x = sd(tbl$FoldChange[10:12])/sqrt(3)
tbl_dotplot = data.frame("Sample" = c("Negative","siRNA_48h"),
                         "Mean" = c(Avg.Negative,Avg.siRNA_48h), 
                         "se" = c(SD.Negative,SD.siRNA_48h))
tbl_dotplot$Sample = factor(c("Negative","siRNA_48h"),
                            levels = c("Negative","siRNA_48h"))
tbl$Group = c(rep(c('Negative',"siRNA_48h"),each = 2))
#tbl.new = tbl[,c(1,14:17)]
#tbl.dotplot = melt(tbl.new)
tbl$Samples = factor(c("Negative_1","Negative_2","siRNA_48h_1","siRNA_48h_2"),
                     levels = c("Negative_1","Negative_2","siRNA_48h_1","siRNA_48h_2"))
#tbl.dotplot$variable = c(rep(c('IronTreat_10mM', 'IronTreat_siRNA', 'IronTreat_Negative', 'IronTreat_Ctrl'), each = 3))

p = ggplot(tbl_dotplot,aes(x = Sample, y = Mean)) + geom_point(aes(size = 4),stat='identity', position = 'dodge', shape = 21)+
  geom_errorbar(aes(ymin = Mean - se, ymax = Mean + se), position =position_dodge(0.9), width = 0.15)+
#  geom_point(aes(x = Group, y = FoldChange), tbl, position = position_dodge(0.9))+
  theme(legend.title = element_blank(), plot.title = element_text(size = 30, hjust = 0.5,margin = margin(0,0,40,0)), 
        panel.background = element_blank(), axis.line = element_line(size = 1, colour = 'black'), 
        legend.position = 'none', axis.text.x = element_text(size = 15, colour = 'black', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15, colour = 'black'), axis.title.y = element_text(size = 20))+ 
  xlab(NULL) +ylab('Fold change\n(IronTreat/Control)')+ 
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.8))+
  labs(title = "Ptp4a3 siRNA 48h (n = 2)")+ scale_fill_grey()
#  geom_signif(y_position=c(1.4), xmin=c(1), xmax=c(2),
#              annotation=c("N.S."), tip_length=0.05, textsize = 5)
  
p

t.test(tbl$FoldChange[1:3],tbl$FoldChange[4:6])
shapiro.test(tbl$FoldChange) #Normality
levene.test(tbl$FoldChange, tbl$Group, location = 'mean') #Homogeniety of Variance Test
var.test(FoldChange ~ Group, data = tbl)
anova_group = aov(FoldChange ~ Group, data = tbl)
bartlett.test(ddCt_FoldChange ~ group, data = group_df)
TukeyHSD(anova_group)
