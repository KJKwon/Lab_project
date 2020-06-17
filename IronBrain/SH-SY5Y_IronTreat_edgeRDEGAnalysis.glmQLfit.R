library(edgeR)
tbl = read.table('SH-SY5Y_IronTreat_count_GeneName.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('SH-SY5Y_IronTreat_rpkm_GeneName.txt', header = TRUE, sep = '\t', row.names = 1)

#Correlation between sample (Standard for undergoing analysis) 
tmp_cor <- dist( 1 - cor(as.matrix(tbl_rpkm), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
plot(tmp_clust, main = 'SH-SY5Y_IronTreat', cex = 1.5, cex.main = 1.5)
#select = apply(tbl_rpkm, 1, function(x) sum(x >= 1) == 9 )
#tbl.clean = tbl[select,]
#tbl_rpkm.clean = tbl_rpkm[select,]

#Ctrl vs 10mM
tbl.Ctrlvs10mM = tbl[,c(1,2,3,7,8,9)]
tbl.Ctrlvs10mM.rpkm = tbl_rpkm[,c(1,2,3,7,8,9)] 
pair3.select = apply(tbl.Ctrlvs10mM.rpkm, 1, function(x) sum(x >= 1) >= 3)
tbl.Ctrlvs10mM.clean = tbl.Ctrlvs10mM[pair3.select,]
Conc.group = factor(c('Ctrl','Ctrl','Ctrl','10mM','10mM','10mM'))
y = DGEList(tbl.Ctrlvs10mM.clean,group = Conc.group)
y = calcNormFactors(y)
plotMDS(y)
#For robust option refer to http://supportupgrade.bioconductor.org/p/79149/#79225
design = model.matrix(~Conc.group)
y = estimateDisp(y, design)
fit = glmQLFit(y, design)
qlf = glmQLFTest(fit, coef = 2)
write.table(topTags(qlf , n = Inf),'SH-SY5Y_IronTreat_Ctrlvs10mM.edgeRQLF_output.txt', sep = '\t',quote = FALSE)
DEG.total = topTags(qlf , n = Inf)$table
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$FDR < 0.05, ]
write.table(DEG.clean,'SH-SY5Y_IronTreat_Ctrlvs10mM.edgeRQLF_output_clean.txt', sep = '\t',quote = FALSE)
#Ctrl vs 2mM vs 10mM (Compare by each concentration)
tbl.Ctrlvs2mMvs10mM = tbl[,c(2,3,4,5,8,9)]
tbl.Ctrlvs2mMvs10mM.rpkm = tbl_rpkm[,c(2,3,4,5,8,9)]
pair2.select = apply(tbl.Ctrlvs2mMvs10mM.rpkm, 1, function(x) sum(x >= 1) >= 2)
tbl.Ctrlvs2mMvs10mM.clean = tbl.Ctrlvs2mMvs10mM[pair2.select,]
Conc.group = factor(c('Ctrl','Ctrl','2mM','2mM','10mM','10mM'))
y = DGEList(tbl.Ctrlvs2mMvs10mM.clean, group = Conc.group)
y = calcNormFactors(y)
design = model.matrix(~0+Conc.group)
y = estimateDisp(y, design)
fit = glmQLFit(y, design)
qlf.Ctrlvs2mM = glmQLFTest(fit, contrast = c(0,1,-1))
qlf.2mMvs10mM = glmQLFTest(fit, contrast = c(1,-1,0))
qlf.Ctrlvs10mM = glmQLFTest(fit, contrast = c(1,0,-1))
write.table(topTags(qlf.Ctrlvs2mM , n = Inf),'SH-SY5Y_IronTreat_Ctrlvs2mM.3pairs.edgeRQLF_output.txt', sep = '\t',quote = FALSE)
DEG.total = topTags(qlf.Ctrlvs2mM , n = Inf)$table
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$FDR < 0.05, ]
write.table(DEG.clean,'SH-SY5Y_IronTreat_Ctrlvs2mM.3pairs.edgeRQLF_output_clean.txt', sep = '\t',quote = FALSE)
write.table(topTags(qlf.2mMvs10mM , n = Inf),'SH-SY5Y_IronTreat_2mMvs10mM.3pairs.edgeRQLF_output.txt', sep = '\t',quote = FALSE)
DEG.total = topTags(qlf.2mMvs10mM , n = Inf)$table
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$FDR < 0.05, ]
write.table(DEG.clean,'SH-SY5Y_IronTreat_2mMvs10mM.3pairs.edgeRQLF_output_clean.txt', sep = '\t',quote = FALSE)
write.table(topTags(qlf.Ctrlvs10mM , n = Inf),'SH-SY5Y_IronTreat_Ctrlvs10mM.3pairs.edgeRQLF_output.txt', sep = '\t',quote = FALSE)
DEG.total = topTags(qlf.Ctrlvs10mM , n = Inf)$table
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$FDR < 0.05, ]
write.table(DEG.clean,'SH-SY5Y_IronTreat_Ctrlvs10mM.3pairs.edgeRQLF_output_clean.txt', sep = '\t',quote = FALSE)
