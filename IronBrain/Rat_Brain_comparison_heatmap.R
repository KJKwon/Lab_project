library(gplots)
library(RColorBrewer)
#Color section
Colors <<- rev(brewer.pal(11,"Spectral"))
Colors <<- colorRampPalette(Colors)(100)
#Remove expression noise
Cleantbl = function(tbl, num){
  select = apply(tbl, 1, function(x) sum(x >= 1) >= num )
  tbl.clean = tbl[select,]
  return(tbl.clean)
}
#Dendrogram (Based on correlation)
Cor.test = function(tbl.clean){
  sample_cor = dist(1 - cor(as.matrix(tbl.clean), method = 'spearman'))
  sample_clust = hclust(sample_cor, method = 'average')
  sample_clust = as.dendrogram(sample_clust)
  plot(sample_clust, main = "Female 6Mo vs 6Wo Left+Right")
}

#log fold change heatmap
ExpNormtoHeatmap = function(tbl.clean, threshold){
  ExpNorm.tbl.clean = tbl.clean
  num.col = length(colnames(ExpNorm.tbl.clean))
  select = apply(ExpNorm.tbl.clean, 1, function(x) sum( x == 0 ) == 0 )
  ExpNorm.tbl.clean = ExpNorm.tbl.clean[select,]
  tmp_logFC = as.vector(apply(ExpNorm.tbl.clean[,(num.col/2 + 1):num.col], 1, mean))
  ExpNorm.tbl.clean.logFC = ExpNorm.tbl.clean[,1:(num.col/2)] / tmp_logFC
  ExpNorm.tbl.clean.logFC = ExpNorm.tbl.clean.logFC[apply(ExpNorm.tbl.clean.logFC, 1, function(x) sum(abs(log(x,base = 2)) < threshold) <= 2),]
  ExpNorm.tbl.clean.logFC = ExpNorm.tbl.clean.logFC[order(apply(ExpNorm.tbl.clean.logFC[,1:(num.col/2)], 1, mean)),]
  k = heatmap.2(as.matrix(log(ExpNorm.tbl.clean.logFC, base = 2)), trace = "none", col = Colors,Rowv = FALSE, Colv = FALSE,
                main = "Female 6Mo vs 6Wo RPKM \nlog2FC heatmap (log2FC > 0.2)")
  return(ExpNorm.tbl.clean.logFC)
}

#Finalize heatmap (Final list selection)
FinalHeatmap = function(tbl.clean, ExpNorm.tbl.clean.logFC){
  ExpNorm.tbl.clean = tbl.clean
  num.col = length(colnames(ExpNorm.tbl.clean))
  ExpNorm.tbl.revise = ExpNorm.tbl.clean[rownames(ExpNorm.tbl.clean) %in% rownames(ExpNorm.tbl.clean.logFC),]
  select = apply(ExpNorm.tbl.revise, 1, function(x) sum( x == 0 ) == 0 )
  ExpNorm.tbl.revise = ExpNorm.tbl.revise[select,]
  ExpNorm.tbl.revise = ExpNorm.tbl.revise[order(apply(ExpNorm.tbl.revise[,1:(num.col/2)], 1, mean)),]
  k = heatmap.2(as.matrix(log(ExpNorm.tbl.revise, base = 2)), trace = "none", col = Colors,scale = "row", Rowv = FALSE, Colv = FALSE, main = "Female 6Mo vs 6Wo heatmap [Age]")
  ExpNorm.tbl.revise.norm = t(as.data.frame(k$carpet))
  #Heatmap filtering
  filter.old = rowSums(ExpNorm.tbl.revise.norm[,1:(num.col/2)] > 0) <= round(num.col*0.1) | rowSums(ExpNorm.tbl.revise.norm[,1:(num.col/2)] < 0) <= round(num.col*0.1)
  filter.young = rowSums(ExpNorm.tbl.revise.norm[,(num.col/2 + 1):num.col] > 0) <= round(num.col*0.1) | rowSums(ExpNorm.tbl.revise.norm[,(num.col/2 + 1):num.col] < 0) <= round(num.col*0.1)
  ExpNorm.tbl.revise.norm = ExpNorm.tbl.revise.norm[filter.old & filter.young,]
  ExpNorm.tbl.revise.norm = ExpNorm.tbl.revise.norm[order(apply(ExpNorm.tbl.revise.norm[,(num.col/2 + 1):num.col], 1, mean)),]
  heatmap.2(as.matrix(ExpNorm.tbl.revise.norm), trace = "none", col = Colors, scale = "none", Rowv = FALSE, Colv = FALSE,
            main = "Male 15Mo vs 6Mo RPKM candidate heatmap \n [log2FC > 0.2 cut]", keysize = 1, lhei = c(1,6), lwid = c(1.2,5),
            margins = c(10,7), cexRow = 0.5)
  
  return(rownames(ExpNorm.tbl.revise.norm))
}

tbl = read.table('SH-SY5Y_IronTreat_Ctrlvs10mM.edgeRQLF_output_clean_RatOrtho.txt', sep = '\t', header = TRUE)
Rat_ortho = unique(tbl$RatOrtho)
Rat_rpkm_tbl = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', sep ='\t', header = TRUE, row.names = 1)

##TPM SECTION START##
#Rat_tpm_tbl = t(t(Rat_rpkm_tbl) / as.vector(apply(Rat_rpkm_tbl, 2, sum))) * (10^6)
#write.table(Rat_tpm_tbl,"Rat_Brain_15M_6M_6W_tpm_GeneName.txt", sep = '\t', quote = FALSE)
#Rat_tpm_tbl = read.table('Rat_Brain_15M_6M_6W_tpm_Left+Right_GeneName.txt', sep = '\t', header = TRUE, row.names = 1)
#Rat_tpm_candi = Rat_tpm_tbl[rownames(Rat_tpm_tbl) %in% Rat_ortho,]
#Rat_tpm_candi.male = Rat_tpm_candi[,1:16]
#Rat_tpm_candi.female = Rat_tpm_candi[,17:36]
#Rat_tpm_male = Rat_tpm_candi[,1:8]
#Rat_tpm_female = Rat_tpm_candi[,9:18]
##TPM SECTION END##

Rat_rpkm_tbl = read.table('Rat_Brain_15M_6M_6W_rpkm_Left+Right_GeneName.txt', sep = '\t', header = TRUE , row.names = 1)
Rat_rpkm_candi = Rat_rpkm_tbl[rownames(Rat_rpkm_tbl) %in% Rat_ortho,]
Rat_rpkm_male = Rat_rpkm_candi[,1:8]
Rat_rpkm_female = Rat_rpkm_candi[,9:18]
Rat_rpkm_male.filtered = Rat_rpkm_male[,3:6]
#Rat_rpkm_candi.male = Rat_rpkm_candi[,1:16]
#Rat_rpkm_candi.female = Rat_rpkm_candi[,17:36]

tmp_clean = Cleantbl(Rat_rpkm_male.filtered,2)
Cor.test(tmp_clean)
tmp_logFC = ExpNormtoHeatmap(tmp_clean, 0.2)
gene.candidate = FinalHeatmap(tmp_clean, tmp_logFC)
length(gene.candidate)
