library(edgeR)
tbl = read.table('SH-SY5Y_IronTreat_24h_MGISeq_count_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('SH-SY5Y_IronTreat_24h_MGISeq_rpkm_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
select = apply(tbl_rpkm[,c(1,2,3,10,11,12)], 1, function(x) sum(x >= 1) >= 3)
tbl_rpkm.clean = tbl_rpkm[select,]
tbl.clean = tbl[select,]
tmp_cor <- dist( 1 - cor(as.matrix(tbl_rpkm.clean), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
par(mar = c(10,4,4,2))
plot(tmp_clust, main = 'SH-SY5Y IronTreat', cex = 1.5, cex.main = 1.5)

#Ctrl vs 10mM
tbl.clean = tbl.clean[,c(1,2,3,10,11,12)]
tbl_rpkm.clean = tbl_rpkm.clean[,c(1,2,3,10,11,12)] 
Conc.group = factor(c('10mM','10mM','10mM','Ctrl','Ctrl','Ctrl'), levels = c('Ctrl','10mM'))
y = DGEList(tbl.clean ,group = Conc.group)
y = calcNormFactors(y)
plotMDS(y)
#For robust option refer to http://supportupgrade.bioconductor.org/p/79149/#79225
design = model.matrix(~Conc.group)
y = estimateDisp(y, design)
fit = glmQLFit(y, design, robust = TRUE)
plotQLDisp(fit)
qlf = glmQLFTest(fit, coef = 2)
write.table(topTags(qlf , n = Inf),'SH-SY5Y_IronTreat_Ctrlvs10mM_MGISeq.edgeRQLF_Robust_geneID_output.txt', sep = '\t',quote = FALSE)
DEG.total = topTags(qlf , n = Inf)$table
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$FDR < 0.05, ]
write.table(DEG.clean,'SH-SY5Y_IronTreat_Ctrlvs10mM.MGISeq.edgeRQLF_Robust_geneID_output_clean.txt', sep = '\t',quote = FALSE)

tbl.Ctrlvs10mM.rpkm.filtered = tbl.Ctrlvs10mM.rpkm[rownames(tbl.Ctrlvs10mM.rpkm) %in% rownames(y$counts), ]
par(mar = c(8,2,2,2))
tmp_cor <- dist( 1 - cor(as.matrix(tbl.clean), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
plot(tmp_clust, main = 'SH-SY5Y_IronTreat_RNAseq_filtered', cex = 1.5, cex.main = 1.5)

#Compare by each concentration
tbl.Ctrlvs2mMvs10mM = tbl.clean[,c(2,3,4,5,7,8)]
tbl.Ctrlvs2mMvs10mM.rpkm = tbl_rpkm.clean[,c(2,3,4,5,7,8)]
pair2.select = apply(tbl.Ctrlvs2mMvs10mM.rpkm, 1, function(x) sum(x >= 1) >= 2)
tbl.Ctrlvs2mMvs10mM.clean = tbl.Ctrlvs2mMvs10mM[pair2.select,]
tbl.Ctrlvs2mMvs10mM.rpkm.clean = tbl.Ctrlvs2mMvs10mM.rpkm[pair2.select,]
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
