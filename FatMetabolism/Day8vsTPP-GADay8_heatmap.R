library(circlize)
library(ComplexHeatmap)
tbl = read.table('Day8vsTPP-GA_Day8.txt', header = TRUE)
GO_pos_reg_fat_diff = c('Adig','Rarres2','Igf1','Cebpa','Sfrp2','Pparg','Noct','Lrp4','Snai2')
GO_reg_tri_biosyn_process = c('Thrsp', 'Dgat2', 'C3', 'Nr1h3', 'Scarb1' ,'Fitm2')
GO_beta_oxy = c('Adipoq','Abcd2', 'Hadhb', 'Ivd', 'Acadm', 'Acad11', 'Echs1', 'Etfdh','Acadvl',
                 'Acads', 'Scp2', 'Crat', 'Acadsb', 'Abcd1', 'Cpt1a', 'Acox2')
#GO_Adipokine = c('Cfd','Retn','Adipoq','Rarres2','Igf1','Tnfrsf21','Igf2','Tnfaip6','Lcn2','Cxcl10','Nampt',
#            'Vegfc','Tnfrsf12a','Ccl5','Ctsl','Ddit3')
#GO_list = GO_Adipokine
GO_neg_reg_to_ER_stress = c('Hspa8','Ppp1r15a','Rack1','Herpud1','Hyou1','Hspa5','Wfs1')
GO_list = c(GO_pos_reg_fat_diff,GO_reg_tri_biosyn_process,GO_beta_oxy,GO_neg_reg_to_ER_stress)
tbl.selected = tbl[match(GO_list,tbl$GeneName),]
rownames(tbl.selected) = tbl.selected$GeneName
tbl.selected = tbl.selected[,c(-1)]
tbl.selected.sc = apply(tbl.selected, 1, function(x){(x-mean(x))/sd(x)})
tbl.selected.sc = t(tbl.selected.sc)
GO_group = c(rep(1,length(GO_pos_reg_fat_diff)),rep(2,length(GO_reg_tri_biosyn_process)),rep(3, length(GO_beta_oxy)),
                   rep(4,length(GO_neg_reg_to_ER_stress)))
col_fun = colorRamp2(c(-max(tbl.selected.sc),0,max(tbl.selected.sc)), c("blue","white","red")) 
Heatmap(as.matrix(tbl.selected.sc), row_order = rownames(tbl.selected), column_order = colnames(tbl.selected), col = col_fun,
        heatmap_width = unit(6, "cm"),heatmap_height = unit(10,"cm"), show_column_names = FALSE, row_split = 
        heatmap_legend_param = list(grid_height = unit(10,"mm")))
