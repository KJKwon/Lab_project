library(edgeR)
tbl = read.table('Rat_Brain_15M_6M_6W_count_GeneName.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', header = TRUE, sep = '\t', row.names = 1)

select <- apply(tbl_rpkm, 1, function(x) sum(x >= 1) >= 16)
tbl_clean = tbl[select,]
tbl_rpkm_clean = tbl_rpkm[select,]
Age.group = factor(c(rep('M15', times = 8), rep('M6', times = 18), rep('W6', times= 10)),
                   levels = c('W6','M6','M15'))
y <- DGEList(tbl, group = Age.group)
y <- calcNormFactors(y)
design <- model.matrix(~0 + Age.group)
colnames(design) <- levels(Age.group)
y <- estimateDisp(y, design)
fit = glmQLFit(y, design, robust = TRUE)
con <- makeContrasts('15Mvs6M' = M15-M6,
                     '6Mvs6W' = M6-W6,
                     '15Mvs6W' = M15-W6, 
                     levels= design)
test <- glmQLFTest(fit, contrast = c(0,-1,1))

topTags(anov)
write.table(topTags(test, n = Inf), 'Rat_Brain_6Mvs6W_count_fastq.trimmed.RAT_ens98_longest_cDNA.BWA.SAMBestCount.glmQLFit.RobustDisp.ANOVA.edgeR.txt',
            sep = '\t', quote = FALSE)
