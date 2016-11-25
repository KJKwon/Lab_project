library(VennDiagram)
grid.newpage()
draw.pairwise.venn(2419,3610,1130,category = c("Human-positive","Mouse-positive"),fill = c("light blue","pink"), lty = rep("blank",2), 
                   alpha = rep(0.5,2), cat.pos= c(180,180), cat.cex = rep(2,2),cex = rep(2,3))
grid.text('(7611)',x = unit(0.22,"npc"),y=unit(0.435,"npc"), draw = TRUE)
grid.text('(6586)',x = unit(0.85,"npc"),y=unit(0.435,"npc"), draw = TRUE)

grid.newpage()
draw.pairwise.venn(2279,3306,810,category = c("Human-negative","Mouse-negative"),fill = c("light blue","pink"), lty = rep("blank",2), 
                   alpha = rep(0.5,2), cat.pos= c(180,180), cat.cex = rep(2,2),cex = rep(2,3))
grid.text('(8024)',x = unit(0.226,"npc"),y=unit(0.435,"npc"), draw = TRUE)
grid.text('(6420)',x = unit(0.834,"npc"),y=unit(0.435,"npc"), draw = TRUE)
