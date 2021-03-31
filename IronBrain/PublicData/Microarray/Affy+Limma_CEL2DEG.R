  library(affy)
library(limma)
library(mouse4302.db)
library(annotate)
cel.path = getwd()
cel.list = list.files(cel.path,full.names = FALSE,pattern = '.CEL')
cel.3wks.list = cel.list[c(1:3,5:7)]
cel.10wks.list = cel.list[c(8:10,12:14)]
eset <- justRMA(filenames = cel.10wks.list)
#colnames(eset) = c("Control_3w_1","Control_3w_2","Control_3w_3",
#                   "Tfr1_null_3w_2","Tfr1_null_3w_3","Tfr1_null_3w_4")
colnames(eset) = c("Tfr1_null_10w_1","Tfr1_null_10w_2","Tfr1_null_10w_3",
                   "Control_10w_2","Control_10w_3","Control_10w_4")
h = hclust(dist(1-cor(exprs(eset))), method = 'complete')
plot(h, main = 'Control vs Tfr1_null 3wks')
ID = featureNames(eset)
Symbol = getSYMBOL(ID, "mouse4302.db")
fData(eset) = data.frame(Symbol = Symbol)
HasSymbol <- !is.na(fData(eset)$Symbol)
eset <- eset[HasSymbol,]

Exp <- factor(c("Tfr1_null","Tfr1_null","Tfr1_null",
                "Control","Control","Control"), levels = c("Control","Tfr1_null"))
design <- model.matrix(~Exp)
fit <- lmFit(eset, design)

fit <- eBayes(fit, trend= TRUE, robust = TRUE)
results <- decideTests(fit)
summary(results)
out = topTable(fit,number = Inf)
write.table(out,'Control_10w_vs_Tfr1_null_10w_GSE66730_filtered_limma.txt', sep = '\t', quote = FALSE)
