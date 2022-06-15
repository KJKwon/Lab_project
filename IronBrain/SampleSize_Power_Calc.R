library(RnaSeqSampleSize)
set.seed(123)
tbl = read.table('GSE66715_Count.csv', sep = '\t', row.names = 1, header = TRUE)
keep = apply(tbl, 1, function(x) sum(x == 0) < 3 )
tbl = tbl[keep,]
ref_data <- est_count_dispersion(tbl, group = c(rep(0,3),rep(1,3)))
est_power_distribution(n=8,f=0.05,rho=2,distributionObject=ref_data,repNumber=100)
tbl_rat <- read.table("Rat_Brain_15M_6M_6W_count_GeneName.txt", sep = '\t', row.names = 1, header =TRUE)
tbl_rat_set <- tbl_rat[,c(1:26)]
our_data <- est_count_dispersion(tbl_rat_set, group = c(rep(0,8),rep(1,18)))
results<- est_power_distribution(n = 16, f = 0.05, rho = 2, distributionObject = our_data, repNumber = 100)
sample_size_distribution(power = 0.8, f = 0.05, distributionObject = our_data, repNumber = 100,
                         showMessage = TRUE, k = 2)
sample_size(power = 0.8, f = 0.05, k = 1, rho = 2, lambda0 = 551.8924 , phi0 = 0.00592)
