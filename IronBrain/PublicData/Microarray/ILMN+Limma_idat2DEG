library(limma)
##Using Non-normalized data
dat <- read.delim("GSE70431_Non-normalized_data.txt", sep = '\t', row.names = 1)
j <- 2*(1:8)
x <- new("EListRaw")
x$E <- as.matrix(dat[,j-1])
x$other$Detection <- as.matrix(dat[,j])
y <- neqc(x)
expressed <- rowSums(y$other$Detection < 0.05) >= 4
y <- y[expressed,]
h = hclust(dist(1-cor(y$E)), method = 'complete')
plot(h)
y$targets <- data.frame('Sample' = c('Mut1','Mut2','Mut3','Mut4','Ctrl1','Ctrl2','Ctrl3','Ctrl4'))
colnames(y) = y$targets$Sample
Exp <- factor(c("Mut","Mut","Mut",
                "Control","Control","Control"), levels = c("Control","Mut"))
design <- model.matrix(~Exp)
fit <- lmFit(y$E[,c(1,3,4,6,7,8)], design)
fit <- eBayes(fit, trend= TRUE, robust = TRUE)
results <- decideTests(fit)
summary(results)
out = topTable(fit,number = Inf)

##Using Raw idat file and bgx
idat.files <- dir(pattern = '.idat')
bgxfile <- "GPL6885_MouseRef-8_V2_0_R0_11278551_A.bgx"
DoubleMut.project <- read.idat(idatfiles = idat.files,bgxfile, annotation = 'Symbol')
DoubleMut.project$other$Detection <- detectionPValues(DoubleMut.project)
DoubleMut.project$targets <- data.frame('Sample' = c('Mut1','Mut2','Mut3','Mut4',
                                                     'Ctrl1','Ctrl2','Ctrl3','Ctrl4'))
DoubleMut.project <- DoubleMut.project[DoubleMut.project$genes$Status == 'regular',]
colnames(DoubleMut.project) <- DoubleMut.project$targets$Sample
expressed <-rowSums(DoubleMut.project$other$Detection < 0.05) >= 4
DoubleMut.project <- DoubleMut.project[expressed,]
h = hclust(dist(1-cor(y$E)), method = 'complete')
plot(h)
