library(edgeR)
mousedata <- read.table("mouse_counted_data_SN_vs_VTA_filt_ensembl.txt",row.names = 1, header = TRUE)
mousedatagroup <- c('SN','SN','SN','SN','SN','SN','VT','VT','VT','VT','VT','VT')
dge <- DGEList(counts = mousedata, group = factor(mousedatagroup))
dge <- calcNormFactors(dge)
design <- read.table("Mouse_sup.txt",row.names = 1, header = TRUE)
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
plotMA(fit, coef = ncol(fit))
adj.p <- p.adjust(fit$p.value[,2], method="fdr")
pos_m_list <- c('Bfsp2', 'A3Galt2', 'Ggt1', 'Catsper3', 'Gsg2', 'Syce1', 'Vsig2', 'Hcrt', 'Bpifb2', 'Mgat4D', 'Steap4')
neg_m_list <- c('Mybl2', 'Galr2', 'Csf2Rb', 'Lrrc23', 'Micalcl', 'Zfp36L1', 'Tap1', 'Mmp28', 'P2Ry2', 'B3Gnt5', 'Nt5E', 'Frmd7', 'Fst', 'Acss3', 'Ripk4', 'Nlrx1', 'Espl1', 'Cabp4', 'Sox21', 'Dpp4', 'Atp7B', 'Copz2', 'Epyc', 'Lgi4', 'Il15', 'Pla2G2F', 'Neurog2', 'Ptch2', 'Cfap161', 'Cmtm3', 'Irf4', 'Pdzd2', 'Notch2', 'Rln3', 'Tnfrsf1A', 'Lad1', 'Cd68', 'Zmynd10', 'Fam83B', 'Hspb2', 'Kctd11', 'Glis3', 'Adgrf2', 'Gpc3', 'Gpc4', 'Rasef', 'Kif12', 'Pdgfrb', 'Plpp4', 'Fyb', 'Cdh24', 'Meox2', 'Cxcl14', 'Zic2', 'Adamtsl4', 'Has2', 'Ackr3', 'Crocc2', 'Lyn', 'Tlr5', 'Col26A1', 'Neurod6', 'Arsj', 'Ranbp2', 'Ets1', 'Mertk', 'Npas4', 'Gpr25', 'Epha1', 'Kdm6B', 'Fgf17', 'Col9A1', 'Kcnj14', 'Prr12', 'Mc3R', 'Sstr4', 'Igsf9', 'Gkn1', 'Rad54B', 'Trpm8', 'Galnt12', 'Esrrb', 'Fam19A3', 'Pof1B', 'Tsku', 'Msc', 'Pdgfb', 'Igdcc3', 'Lcorl', 'Celsr1', 'Pax7', 'Dsg2', 'Fbxo40', 'Lrrk1', 'Itpr2', 'Ncor2', 'Hivep1', 'Hs3St3B1', 'Grin3B', 'Myocd', 'Atp6V1C2', 'Chsy3', 'Col12A1', 'Mst1R', 'Fam150B', 'Pxdc1', 'Lrrn4', 'Hnf1B', 'Cdc42Bpg', 'Hsd17B14', 'Pou2F2', 'Casp6', 'Unc13C', 'Pdlim1', 'Hes5', 'Smo', 'Lpl', 'Galnt3', 'Cd247', 'Vipr1', 'Neurod2', 'Neurod4', 'Nrap', 'Cpa4', 'Dach2', 'Irx1', 'Muc6', 'Scn11A', 'Trim15', 'Npffr2', 'Masp1', 'Hspg2', 'Tmem204', 'Sox1', 'Eya1', 'Fxyd3', 'Nek10', 'Ptx3', 'Klf2', 'Snai2', 'Cnih3', 'Meig1', 'Elovl3', 'Heyl', 'Dusp5', 'Rbl1', 'Fam83F', 'Tstd1', 'Vwa5B1', 'Dnah2', 'Eva1A', 'Adamts8', 'Prokr2', 'Map3K14', 'Cenpa', 'Pabpc4L', 'Fgd5', 'Mc4R', 'Aim1', 'Nr4A1', 'Dchs1', 'Slc14A2', 'Yap1', 'Eda2R', 'Cd4', 'Edn1', 'Enpp3', 'Dkkl1')
neg_m_point <- c()
G <- rownames(mousedata)
for( i in 1:length(adj.p)){
  for(j in 1:length(pos_m_list)){
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
output <-topTable(fit, number = Inf, coef=ncol(design))
write.table(output, "Mouse_sig_output_ensembl.txt",sep = '\t')
