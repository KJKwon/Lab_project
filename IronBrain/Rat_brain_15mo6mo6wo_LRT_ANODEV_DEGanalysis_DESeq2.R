library(DESeq2)
tbl = read.table('Rat_Brain_15M_6M_6W_count_GeneName.txt', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm_GeneName.txt', header = TRUE, sep = '\t', row.names = 1)
select <- apply(tbl_rpkm, 1, function(x) sum(x >= 1) >= 8)
tbl_clean = tbl[select,]
tbl_rpkm_clean = tbl_rpkm[select,]
coldata <- data.frame(Age = factor(c(rep('M15',times = 8), rep('M6',times = 18),rep('W6', times = 10))))
rownames(coldata) = colnames(tbl)
coldata$Age = factor(coldata$Age, levels = c('M6','W6','M15'))
dds <- DESeqDataSetFromMatrix(countData = tbl_clean, colData = coldata, design = ~Age)
#sigLRT_vst <- vst(dds)
#sigLRT_vst_mat <- assay(sigLRT_vst)
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res_15Mvs6M <- results(dds, name = 'Age_M15_vs_M6')
res_6Mvs6W <- results(dds, name = 'Age_M6_vs_W6')
write.table(res_15Mvs6M, 'Rat_Brain_15M_6M_6W.fastq.trimmed.RAT_ens98_longest_cDNA.BWA.DESeq2_LRT_15mo_vs_6mo_GeneName.txt',
            quote = FALSE, sep = '\t')
