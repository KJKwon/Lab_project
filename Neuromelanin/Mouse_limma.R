library(edgeR)
library(limma)
comp = read.table('Mouse_count_data.txt',row.names = 1, header = TRUE,sep = '\t')
rownames(comp) <- comp$Names
comp <- subset(comp, select = -c(Names))
comp <- comp[rowSums(comp)>1,]
dge <- DGEList(counts = comp)
dge <- calcNormFactors(dge)
sup = read.table('Mouse_sup.txt',row.names = 1, header = TRUE, sep = '\t')
type = factor(sup$Types, levels = c("VT", "SN"))
design <- model.matrix(~type)
v <- voom(dge, design)
fit <- lmFit(v,design)
fit <- eBayes(fit)
plotMA(fit, coef = ncol(fit))
adj.p <- p.adjust(fit$p.value[,2], method="fdr")
o <- c()
for (i in 1:length(adj.p)){
  if (adj.p[i] < 0.05){
    if (fit$Amean[i] > 0){
    o <- c(o,i) 
    }
  }
}
o
x <- fit$Amean
y <- fit$coefficients[,2]
smoothScatter(x,y)
G <- rownames(mousedata)
text(x[o],y[o],labels =G[o], cex = 0.5)
top <- topTable(fit,coef =1, number = Inf)
write.table(top,file="Mouse_sig_output_false.txt",sep = "\t")
