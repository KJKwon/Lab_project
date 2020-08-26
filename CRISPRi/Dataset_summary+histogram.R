library(ggplot2)
tbl = read.table('ars1-D_vs_2-D_dataset.txt', header = TRUE, sep = '\t')
tbl.clean = tbl[,c(1,3,5)]
colnames(tbl.clean) = c('GeneName','ars1_D','ars2_D')
rep1 = tbl.clean[,2]
rep2 = tbl.clean[,3]
rep1.IQR = IQR(rep1)
rep2.IQR = IQR(rep2)
rep1.range = c(quantile(rep1)[[2]]-1.5*rep1.IQR, quantile(rep1)[[4]]+1.5*rep1.IQR)
rep2.range = c(quantile(rep2)[[2]]-1.5*rep2.IQR, quantile(rep2)[[4]]+1.5*rep2.IQR)
plot.ready = data.frame(value = c(tbl.clean$ars1_D, tbl.clean$ars2_D),
                        dataset = c(rep(colnames(tbl.clean)[2:3], each = nrow(tbl.clean))))
plot.line = data.frame(IQR = c(rep1.range, rep2.range), datatype = c(rep(c('ars1_D.IQR','ars2_D.IQR'), each = 2)))
ggplot(plot.ready, aes(x = value)) + 
  geom_histogram(aes(fill = dataset),bins = 100, position = 'identity', alpha = 0.5)+
  geom_vline(data = plot.line, aes(xintercept = IQR, color = datatype), linetype = 'dashed')+
  theme_classic() + theme(legend.title = element_blank(),legend.position = "bottom",plot.title = element_text(hjust=0.5),
                          text = element_text(size = 15)) +
  ggtitle("Ars, D") + scale_y_continuous(expand = c(0,0), limits = c(0,1800))

rep.inter = intersect(as.vector(tbl.clean[tbl.clean$ars1_D < rep1.range[1],]$GeneName),as.vector(tbl.clean[tbl.clean$ars2_D < rep2.range[1],]$GeneName)) 
rep.uni = union(as.vector(tbl.clean[tbl.clean$ars1_D < rep1.range[1],]$GeneName),as.vector(tbl.clean[tbl.clean$ars2_D < rep2.range[1],]$GeneName))
Group_D.inter = data.frame(GeneName = rep.inter, GroupD = c(rep('Y',length(rep.inter))))
Group_D.uni = data.frame(GeneName = rep.uni, GroupD = c(rep('Y',length(rep.uni))))

#Total
tbl_total.inter = Reduce(function(x,y) merge(x,y, by = 'GeneName',all = TRUE), list(Group_A.inter, Group_B.inter, Group_C.inter, Group_D.inter, Group_Ref.inter))
tbl_total.inter = as.data.frame(tbl_total.inter, stringsAsFactors = FALSE)
tbl_total.inter = as.matrix(tbl_total.inter)
tbl_total.inter[is.character(tbl_total.inter) & is.na(tbl_total.inter)] = 'N'
tbl_total.inter = as.data.frame(tbl_total.inter)
Outlier_count = apply(tbl_total.inter, 1, function(x) sum(x == 'Y'))
tbl_total.inter = cbind(tbl_total.inter, Outlier_count)
tbl_total.inter = tbl_total.inter[order(tbl_total.inter$GeneName),]
write.table(tbl_total.inter, 'Ars_outlier.intersection.txt', quote = FALSE, sep ='\t', row.names = FALSE)

tbl_total.uni = Reduce(function(x,y) merge(x,y, by = 'GeneName',all = TRUE), list(Group_A.uni, Group_B.uni, Group_C.uni, Group_D.uni, Group_Ref.uni))
tbl_total.uni = as.data.frame(tbl_total.uni, stringsAsFactors = FALSE)
tbl_total.uni = as.matrix(tbl_total.uni)
tbl_total.uni[is.character(tbl_total.uni) & is.na(tbl_total.uni)] = 'N'
tbl_total.uni = as.data.frame(tbl_total.uni)
Outlier_count = apply(tbl_total.uni, 1, function(x) sum(x == 'Y'))
tbl_total.uni = cbind(tbl_total.uni, Outlier_count)
tbl_total.uni = tbl_total.uni[order(tbl_total.uni$GeneName),]
write.table(tbl_total.uni, 'Ars_outlier.union.txt', quote = FALSE, sep ='\t', row.names = FALSE)

