tbl<-read.table('Huebner201609_XENLAtxV2.best_indiv_rpkm.no_batch_humanized_filtered.csv',header=TRUE,sep='\t',row.names=1)
tbl.human <- tbl[-c(1)]
tbl_hgene <- t(tbl.human)
tbl_hgene.corr <- 1-cor(tbl_hgene)
h <- hclust(as.dist(tbl_hgene.corr))
group <- cutree(h,7)
group_cut<-function(group,group_subset){
  group_list <- group[group_subset == TRUE]
  return(group_list)
}
for(x in 1:7){
  group_sub = group == x
  g_sub = group_cut(group,group_sub)
}
