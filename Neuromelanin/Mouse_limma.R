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
#plotMA(fit, coef = ncol(fit))
adj.p <- p.adjust(fit$p.value[,2], method="fdr")
pos_m_list <- c('BFSP2', 'A3GALT2', 'GGT1', 'CATSPER3', 'GSG2', 'SYCE1', 'VSIG2', 'HCRT', 'BPIFB2', 'MGAT4D', 'STEAP4')
neg_m_list <- c('MYBL2', 'GALR2', 'CSF2RB', 'LRRC23', 'MICALCL', 'ZFP36L1', 'TAP1', 'MMP28', 'P2RY2', 'B3GNT5', 'NT5E', 'FRMD7', 'FST', 'ACSS3', 'RIPK4', 'NLRX1', 'ESPL1', 'CABP4', 'SOX21', 'DPP4', 'ATP7B', 'COPZ2', 'EPYC', 'LGI4', 'IL15', 'PLA2G2F', 'NEUROG2', 'PTCH2', 'CFAP161', 'CMTM3', 'IRF4', 'PDZD2', 'NOTCH2', 'RLN3', 'TNFRSF1A', 'LAD1', 'CD68', 'ZMYND10', 'FAM83B', 'HSPB2', 'KCTD11', 'GLIS3', 'ADGRF2', 'GPC3', 'GPC4', 'RASEF', 'KIF12', 'PDGFRB', 'PLPP4', 'FYB', 'CDH24', 'MEOX2', 'CXCL14', 'ZIC2', 'ADAMTSL4', 'HAS2', 'ACKR3', 'CROCC2', 'LYN', 'TLR5', 'COL26A1', 'NEUROD6', 'ARSJ', 'RANBP2', 'ETS1', 'MERTK', 'NPAS4', 'GPR25', 'EPHA1', 'KDM6B', 'FGF17', 'COL9A1', 'KCNJ14', 'PRR12', 'MC3R', 'SSTR4', 'IGSF9', 'GKN1', 'RAD54B', 'TRPM8', 'GALNT12', 'ESRRB', 'FAM19A3', 'POF1B', 'TSKU', 'MSC', 'PDGFB', 'IGDCC3', 'LCORL', 'CELSR1', 'PAX7', 'DSG2', 'FBXO40', 'LRRK1', 'ITPR2', 'NCOR2', 'HIVEP1', 'HS3ST3B1', 'GRIN3B', 'MYOCD', 'ATP6V1C2', 'CHSY3', 'COL12A1', 'MST1R', 'FAM150B', 'PXDC1', 'LRRN4', 'HNF1B', 'CDC42BPG', 'HSD17B14', 'POU2F2', 'CASP6', 'UNC13C', 'PDLIM1', 'HES5', 'SMO', 'LPL', 'GALNT3', 'CD247', 'VIPR1', 'NEUROD2', 'NEUROD4', 'NRAP', 'CPA4', 'DACH2', 'IRX1', 'MUC6', 'SCN11A', 'TRIM15', 'NPFFR2', 'MASP1', 'HSPG2', 'TMEM204', 'SOX1', 'EYA1', 'FXYD3', 'NEK10', 'PTX3', 'KLF2', 'SNAI2', 'CNIH3', 'MEIG1', 'ELOVL3', 'HEYL', 'DUSP5', 'RBL1', 'FAM83F', 'TSTD1', 'VWA5B1', 'DNAH2', 'EVA1A', 'ADAMTS8', 'PROKR2', 'MAP3K14', 'CENPA', 'PABPC4L', 'FGD5', 'MC4R', 'AIM1', 'NR4A1', 'DCHS1', 'SLC14A2', 'YAP1', 'EDA2R', 'CD4', 'EDN1', 'ENPP3', 'DKKL1')
pos_m_point <- c()
neg_m_point <- c()
G <- rownames(mousedata)
for( i in 1:length(adj.p)){
  for(j in 1:length(pos_h_list)){
    p_m_name = pos_m_list[j]
    if(p_m_name == G[i]){
      pos_m_point <-c(pos_m_point, i)
    }
  }
  for(k in 1:length(neg_m_list)){
    n_m_name = neg_m_list[k]
    if(n_m_name == G[i]){
      neg_m_point <-c(neg_m_point,i)
    }
  }
}
x <- fit$Amean
y <- fit$coefficients[,2]
smoothScatter(x,y,xlab="Average log expression", ylab = "log fold change", main ="")
points(x[pos_m_point],y[pos_m_point],pch = 19)
points(x[neg_m_point],y[neg_m_point],pch = 19, col = "red")
text(x[o],y[o],labels =G[o], cex = 0.5)
top <- topTable(fit,coef =1, number = Inf)
write.table(top,file="Mouse_sig_output_false.txt",sep = "\t")
