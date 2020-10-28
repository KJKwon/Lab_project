library(edgeR)
tbl = read.table('Rat_Brain_15M_6M_6W_count_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneID.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_male = tbl[,1:16]
tbl_female = tbl[,17:36]
tbl_rpkm_male = tbl_rpkm[,1:16]
tbl_rpkm_female =tbl_rpkm[,17:36]

#Male case
select = apply(tbl_rpkm_male, 1, function(x) sum(x >= 1) >= 8 )
tbl_male.clean = tbl_male[select,]
tbl_rpkm_male.clean = tbl_rpkm_male[select,]
tmp_cor <- dist( 1 - cor(as.matrix(tbl_rpkm_male), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
plot(tmp_clust, cex = 1.5, cex.main = 1.5)
Age.group = factor(c(rep('15M', times = 8), rep('6M', times = 8)), levels = c('6M','15M'))
#y = DGEList(tbl_male.clean, group = Age.group)
y = DGEList(tbl_male, group = Age.group)
y = calcNormFactors(y)
plotMDS(y) 
#For robust option refer to http://supportupgrade.bioconductor.org/p/79149/#79225
design = model.matrix(~Age.group)
y = estimateDisp(y, design)
fit = glmQLFit(y, design, robust = TRUE)
qlf = glmQLFTest(fit, coef = 2)
write.table(topTags(qlf , n = Inf),'15Month_vs_6Month_selected.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.glmQLFit.RobustDisp.edgeR.txt', sep = '\t',quote = FALSE)
DEG.total = topTags(qlf , n = Inf)$table
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$FDR < 0.05, ]
write.table(DEG.clean,'15Month_vs_6Month_selected.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.glmQLFit.RobustDisp.edgeR_clean.txt', sep = '\t',quote = FALSE)
tbl.male.filtered.rpkm = tbl_rpkm_male[rownames(tbl_rpkm_male) %in% rownames(y$counts), ]
tmp_cor <- dist( 1 - cor(as.matrix(tbl.male.filtered.rpkm), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
plot(tmp_clust, main = 'SH-SY5Y_IronTreat_filtered', cex = 1.5, cex.main = 1.5)

#Female case 
##Step 1) Clean Up data
select = apply(tbl_rpkm_female, 1, function(x) sum(x >= 1) >= 10 )
tbl_female.clean = tbl_female[select,]
tbl_rpkm_female.clean = tbl_rpkm_female[select,]
##Step 2) Correlation check
tmp_cor <- dist( 1 - cor(as.matrix(tbl_rpkm_female.clean), method = 'spearman'))
tmp_clust <- hclust(tmp_cor, method="average")
tmp_clust = as.dendrogram(tmp_clust)
plot(tmp_clust, cex = 1.5, cex.main = 1.5, main = "Female Rat Substantia nigra 6Mo vs 6Wo")
##Step 3) DEG analysis
Age.group = factor(c(rep('6M', times = 10), rep('6W', times = 10)), levels = c('6W','6M'))
y = DGEList(tbl_female.clean, group = Age.group)
y = calcNormFactors(y)
plotMDS(y) 
#For robust option refer to http://supportupgrade.bioconductor.org/p/79149/#79225
design = model.matrix(~Age.group)
y = estimateDisp(y, design)
plotBCV(y)
fit = glmQLFit(y, design, robust = TRUE)
plotQLDisp(fit)
qlf = glmQLFTest(fit, coef = 2)
write.table(topTags(qlf , n = Inf),'6Month_vs_6Week.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.glmQLFit.Robust.edgeR.txt', sep = '\t',quote = FALSE)
DEG.total = topTags(qlf , n = Inf)$table
DEG.clean = DEG.total[abs(DEG.total$logFC) > 1 & DEG.total$FDR < 0.05, ]
write.table(DEG.clean,'6Month_vs_6Week.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.glmQLFit.Robust.edgeR_geneID_clean.txt', sep = '\t',quote = FALSE)
