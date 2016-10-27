tbl <- read.table('Huebner201609_XENLAtxV2.best_indiv_rpkm_batch_removed.txt',header = TRUE, row.names = 1)
group <- rep(1:10,2)
tbl_col = colnames(tbl)
tbl_row = c()
out_name = 'Huebner201609_XENLAtxV2.best_indiv_rpkm_batch_removed_filterered.txt'
for(i in 1:nrow(tbl)){
  B <- c(tbl[i,2],tbl[i,3],tbl[i,4],tbl[i,5],tbl[i,6],tbl[i,7],tbl[i,8],tbl[i,9],tbl[i,10],tbl[i,11])
  C <- c(tbl[i,14],tbl[i,15],tbl[i,16],tbl[i,17],tbl[i,18],tbl[i,19],tbl[i,20],tbl[i,21],tbl[i,22],tbl[i,23])
  BC <- c(B,C)
  break
  group_df <- data.frame(BC,group)
  group_df <- transform(group_df, group = factor(group))
  out <- summary(aov(formula = BC~group, data = group_df))
  p_val = out[[1]]$`Pr(>F)`[[1]]
  if(p_val < 0.05){
    tbl_row[i] <- rownames(tbl[1,])
  }
}
