library(circlize)
library(ComplexHeatmap)
tbl = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneName.txt', header = TRUE)
Iron_related_gene.tbl = read.table('Iron_related_GeneList.txt')
gene_name = Iron_related_gene.tbl$V1
alias_name = Iron_related_gene.tbl$V2
tbl.selected = tbl[match(gene_name,tbl$geneID),]
rownames(tbl.selected) = alias_name
tbl.selected = tbl.selected[,c(-1)]
tbl.selected.sc = apply(tbl.selected, 1, scale)
tbl.selected.sc = as.data.frame(t(tbl.selected.sc))
tbl.selected.sc = tbl.selected.sc[,c(ncol(tbl.selected.sc):1)]
colnames(tbl.selected.sc) = c('Control 1', 'Control 2', 'Control 3', 'FeCl2 500uM 1', 'FeCl2 500uM 2', 'FeCl2 500uM 3',
                              'FeCl2 2mM 1', 'FeCl2 2mM 2','FeCl2 2mM 3', 'FeCl2 10mM 1', 'FeCl2 10mM 2', 'FeCl2 10mM 3')
col_fun = colorRamp2(c(-3,0,3), c("blue","white","red"))
tbl.selected.ordered = tbl.selected.sc[order(apply(tbl.selected.sc[,c(10,11,12)],1,mean), decreasing = TRUE),]
lgd_zscore = Legend(title = "Row z-score", col = col_fun)
ha = HeatmapAnnotation("IronConc." = rep(c(2,3,4,5), each = 3), show_legend = FALSE, 
                       gp = gpar(col = "black"))
ht = Heatmap(as.matrix(tbl.selected.ordered), row_order = rownames(tbl.selected.ordered), column_order = colnames(tbl.selected.ordered),
        col = col_fun,heatmap_width = unit(8, "cm"),heatmap_height = unit(12,"cm"), show_column_names = TRUE, 
        heatmap_legend_param = list(title = "Row z-score", title_position= "leftcenter-rot" ,grid_height = unit(10,"mm")),
        top_annotation = ha)
draw(ht)
