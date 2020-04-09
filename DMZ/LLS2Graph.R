temp.tbl = read.table('Huebner201609_XENLAtxV2_DMZ.indiv_p_rpkm_cor.clean.90LSvalue.txt', sep = '\t',header = TRUE)
temp.tbl.clean = temp.tbl[temp.tbl$Mean.corr > 0.95,]
plot(temp.tbl.clean$Mean.corr, temp.tbl.clean$Likelihood.Score, main = 'DMZ Likelihood score', pch = 20, ylim = c(-3,3),
     xlab = "Mean correlation", ylab = "LLS")
#points(temp.tbl$Mean.corr, temp.tbl$Random.Score, pch = 20, col = "gray")
loess_train = loess(Likelihood.Score ~ Mean.corr, temp.tbl.clean, span = 0.25, control = loess.control(surface = "direct"))
out = predict(loess_train, temp.tbl.clean$Mean.corr)
lines(temp.tbl.clean$Mean.corr, out, col = "red", lwd = 2)
loess_train_rand = loess(Random.Score ~ Mean.corr, temp.tbl.clean, span = 0.25, control = loess.control(surface = "direct"))
out_rand = predict(loess_train_rand, temp.tbl.clean$temp.tbl.clean)
lines(temp.tbl.clean$Mean.corr, out_rand, col = "blue", lwd =2, lty = 2 )
legend(0.95, 4, legend = c('Data regression', 'Random regression'), col = c("red", "blue"), lty = 1:2, cex = 1)
