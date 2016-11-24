library(VennDiagram)
grid.newpage()
draw.pairwise.venn(2419,3610,748,category = c("Human-negative","Mouse-negative"),fill = c("light blue","pink"), lty = rep("blank",2), 
                   alpha = rep(0.5,2), cat.pos= c(180,180), cat.cex = rep(2.5,2),cex = rep(2,3))
