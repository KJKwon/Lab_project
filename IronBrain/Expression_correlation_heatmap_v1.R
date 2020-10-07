library(ComplexHeatmap)
library(RColorBrewer)
tbl.rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', row.names = 1, header = TRUE)
male.filter = apply(tbl.rpkm[,1:16], 1, function(x) sum(x >= 1) >= 8)
female.filter = apply(tbl.rpkm[,17:ncol(tbl.rpkm)], 1, function(x) sum(x >= 1) >= 10)
total.filter = male.filter & female.filter
tbl.rpkm.clean = tbl.rpkm[total.filter,]
colnames(tbl.rpkm.clean) = c(gsub('^Rat_15M','M_15Mo',colnames(tbl.rpkm.clean)[1:8]),gsub('^Rat_6M','M_6Mo',colnames(tbl.rpkm.clean)[9:16]),
                             gsub('^Rat_6M','F_6Mo',colnames(tbl.rpkm.clean)[17:26]),gsub('^Rat_6W','F_6Wo',colnames(tbl.rpkm.clean)[27:36]))
tbl.rpkm.cor = cor(tbl.rpkm.clean, method = 'spearman')
sample.name = colnames(tbl.rpkm.clean)
sample.group = rep(c('Male_15Mo','Male_6Mo','Female_6Mo','Female_6Wo'), c(8,8,10,10))
#tbl.sample = data.frame(sample.name,sample.group)
#sample.group = sample.group[row_order(ht.main)]
bar.group = HeatmapAnnotation(group = sample.group, annotation_name_side = 'left',
                              col = list(group = c("Male_15Mo" = "#F46D43", "Male_6Mo" = "#E6F598",
                                                   "Female_6Mo" = "#ABDDA4", "Female_6Wo" = "#3288BD")),
                              annotation_legend_param = list(grid_height = unit(0.5, "cm"), title = "", grid_with = unit(0.5, "cm")))
ht.main = Heatmap(as.matrix(tbl.rpkm.cor), clustering_distance_rows = "spearman", top_annotation = bar.group,
        clustering_distance_columns = "spearman", 
        heatmap_legend_param = list(legend_height = unit(5, "cm"),title = "Spearman correlation", title_position = "leftcenter-rot",
                                                                              grid_width = unit(0.5,"cm")))
draw(ht.main, merge_legend = TRUE)
