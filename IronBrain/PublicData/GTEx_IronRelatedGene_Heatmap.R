library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
tbl_tpm = read.table('GTExV8_SN_gene_tpm.gct', header = TRUE, check.names = FALSE)
tbl_tpm = tbl_tpm[,c(-1)]
#f_candi = read.table('SH-SY5Y_IronTreat_Controlvs10mM.MGISeq.edgeRQLF_Robust_GeneName_output_clean.txt', header = TRUE, row.names = 1)
#gene_candi = rownames(f_candi)
gene_candi = c('TFRC','SLC40A1','SLC11A2','ACO1','IREB2','FTL','FTH1')
#gene_candi = c('GADD45B','PTP4A3','SULT1A1','VEGFB')
tbl_tpm.selected = tbl_tpm[tbl_tpm$Description %in% gene_candi,]
rownames(tbl_tpm.selected) = tbl_tpm.selected$Description
tbl_tpm.selected = t(tbl_tpm.selected[,c(-1)])
tbl_pheno = read.table('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.full.txt', header = TRUE,row.names = 1)
tbl.complete = merge(tbl_tpm.selected,tbl_pheno, by = intersect(0,0))
#tbl.complete = tbl.complete[tbl.complete$SEX == 'M',]
tbl.complete  = tbl.complete[,c(-1)]
tbl.temp = apply(tbl.complete[,c(1:length(gene_candi))], 2, scale)
tbl.temp = as.data.frame(tbl.temp)
tbl.temp = cbind(tbl.temp, AGE = tbl.complete$AGE)
new.tbl = data.frame()
for(tmp.age in sort(unique(tbl.temp$AGE))){
  tmp.tbl = tbl.temp[tbl.temp$AGE == tmp.age,]
  sum.order = order(tmp.tbl[,1])
  for(i in 2:(length(colnames(tmp.tbl))-1)){
    sum.order = sum.order + order(tmp.tbl[,i])
  }
  tmp.tbl = tmp.tbl[order(sum.order, decreasing = TRUE),]
  if(length(new.tbl) == 0){
    new.tbl = tmp.tbl
  }
  else{
    new.tbl = rbind(new.tbl,tmp.tbl)
  }
}
Age.group = new.tbl$AGE
var = brewer.pal(6, name = "Greys")
col.group = c(rep(var[1], length(Age.group[Age.group == '20-29'])),
              rep(var[2], length(Age.group[Age.group == '30-39'])),
              rep(var[3], length(Age.group[Age.group == '40-49'])),
              rep(var[4], length(Age.group[Age.group == '50-59'])),
              rep(var[5], length(Age.group[Age.group == '60-69'])),
              rep(var[6], length(Age.group[Age.group == '70-79'])))
col = setNames(var,levels(Age.group))
bar.group = HeatmapAnnotation(Age = Age.group, col = list(Age = col),
                              annotation_legend_param = list(grid_height = unit(0.5, "cm"), 
                                                             title = "AGE", grid_width = unit(0.5, "cm"),
                                                             row_names_gp = gpar(fontsize = 12)))
tbl.mat = t(as.matrix(new.tbl[,1:length(gene_candi)]))
dend = cluster_within_group(tbl.mat, Age.group)
p = Heatmap(tbl.mat, column_dend_reorder = FALSE, show_column_dend = FALSE, show_column_names = FALSE,
            cluster_columns = FALSE,top_annotation = bar.group, row_names_gp = gpar(fontsize = 15),
            heatmap_legend_param = list(legend_height = unit(3, "cm"),title = "Row z-score",
                                        grid_width = unit(0.5,"cm"),title_gp = gpar(fontsize = 14),
                                        legend_width = unit(1,"cm"),labels_gp = gpar(fontsize = 13),
                                        title_position = 'leftcenter-rot'))
draw(p)
  #top_annotation = bar.group
h = column_dend(p)
