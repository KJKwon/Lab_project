library(edgeR)
tbl = read.table('Rat_Brain_15M_6M_6W_count.csv', header = TRUE, sep = '\t', row.names = 1)
Outlier = tbl[,23, drop = FALSE]
tbl = tbl[,-c(23)]
Sample = colnames(tbl)
Group = c('E','A','D','D','D','A','D','A','A','E','E','E','C','C','A','A','A',rep('C', times = 8), rep('B', times = 10))
Age = c(rep('15M', times = 8), rep('6M', times = 17), rep('6W', times = 10))
target = data.frame(Sample, Group, Age)
