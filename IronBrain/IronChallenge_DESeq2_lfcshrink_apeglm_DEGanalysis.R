library(DESeq2)
library(DEGreport)
tbl = read.table('SH-SY5Y_IronChallenge_4_lineages.human_ens98_longest_bwa_mem_count_GeneName.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('SH-SY5Y_IronChallenge_4_lineages.human_ens98_longest_bwa_mem_rpkm_GeneName.txt', header = TRUE, sep = '\t', row.names = 1)
coldata <- data.frame(Age = factor(c(rep('Control', times = 3), rep('1mM', times = 3), 
                                     rep('2mM', times= 3),rep('3mM', times= 3))))
rownames(coldata) = colnames(tbl)
coldata$Age = factor(coldata$Age, levels = c('Control','2mM','1mM','3mM'))
dds <- DESeqDataSetFromMatrix(countData = tbl, colData = coldata, design = ~Age)
keep <- rowSums(counts(dds)) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)
resLFC.M15vsM6 <- lfcShrink(dds, coef = "Age_3mM_vs_Control", type = "apeglm")
write.table(resLFC.M15vsM6, 'SH-SY5Y_IronChallenge_3mMvsCtrl_NextSeq_DESeq2_lfcShrink_apeglm_GeneName_output.txt',
            quote = FALSE, sep = '\t')
