library(edgeR)
tbl = read.table('Rat_Brain_15M_6M_6W_count.csv', header = TRUE, sep = '\t', row.names = 1)
tbl_rpkm = read.table('Rat_Brain_15M_6M_6W_rpkm.csv', header = TRUE, sep = '\t', row.names = 1)
select = apply(tbl_rpkm, 1, function(x) sum(x >= 1) == 36 )
Outlier = tbl[,23, drop = FALSE]
tbl = tbl[,-c(23)]
Sample = colnames(tbl)
Group = c('E','A','D','D','D','A','D','A','A','E','E','E','C','C','A','A','A',rep('C', times = 8), rep('B', times = 10))
Age = c(rep(c('15M','6M','6W'),c(8,17,10)))
target = data.frame(Sample, Group, Age)
group = factor(paste0(target$Group, ".", target$Age))
tbl.filtered = tbl[select,]
y = DGEList(tbl.filtered,group = group)
y = calcNormFactors(y)

# Grouping plot section
colors = c('red','red','blue','green','blueviolet','gray0','gray0')
points = c(0,1,2,1,0,0,1)
plotMDS(y, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)

#Validation with dendrogram
tmp_cor = dist(1-cor(cpm(y$counts), method = 'spearman'))
tmp_clust = hclust(tmp_cor, method = 'average')
tmp_plot = as.dendrogram(tmp_clust)
plot(tmp_plot)

#DGE analysis
design = model.matrix(~0 + group)
colnames(design) = levels(group)
y = estimateDisp(y, design, robust = TRUE)
plotBCV(y)
fit = glmQLFit(y, design, robust = TRUE)
con = makeContrasts(C.6M - B.6W, levels = design)
qlf = glmQLFTest(fit, contrast = con)
write.table(topTags(qlf , n = Inf),'groupC_6M_up_groupB_6W_down.edgeRQLF_output.txt', sep = '\t',quote = FALSE)
