library(edgeR)
tbl <- read.table('Huebner201609_XENLAtx_table_B_clean.txt',sep = '\t', head = TRUE,row.names = 1)
DMZ_tbl <- tbl[2:11]
DMZ_compare <- DMZ_tbl[c(3,4,9,10)]
DMZ_comp_group <- c('Early','Early','Late','Late')

new_tbl <- DGEList(counts = DMZ_compare, group = DMZ_comp_group)
new_tbl <- calcNormFactors(new_tbl)
new_tbl <- estimateCommonDisp(new_tbl)
new_tbl <- estimateTagwiseDisp(new_tbl)


write.table(topTags(exactTest(new_tbl,pair = c('Early','Late')),n=Inf),"Hubner201609_XENLAtx_B_Sig_test.txt",sep='\t',quote = FALSE)
write.table(topTags(exactTest(new_tbl,pair = c('Early','Late')),n=Inf),"Hubner201609_XENLAtx_C_Sig_test.txt",sep='\t',quote = FALSE)

B_output <- topTags(exactTest(new_tbl,pair = c('Early','Late')),n=Inf)
C_output <- topTags(exactTest(new_tbl,pair = c('Early','Late')),n=Inf)
B_sig_pos <- subset(B_output$table ,FDR < 0.05 & logFC > 1 )[,c(1,4)]
C_sig_pos <- subset(C_output$table, FDR < 0.05 & logFC > 1 )[,c(1,4)]
B_sig_neg <- subset(B_output$table ,FDR < 0.05 & logFC < -1 )[,c(1,4)]
C_sig_neg <- subset(C_output$table, FDR < 0.05 & logFC < -1 )[,c(1,4)]

total_gene_pos <- intersect(rownames(B_sig_pos),rownames(C_sig_pos))
total_gene_neg <- intersect(rownames(B_sig_neg),rownames(C_sig_neg))

B_sig_pos <- B_sig_pos[c(total_gene_pos),]
C_sig_pos <- C_sig_pos[c(total_gene_pos),]
B_sig_neg <- B_sig_neg[c(total_gene_neg),]
C_sig_neg <- C_sig_neg[c(total_gene_neg),]

#Aim : make final table with logFC,FDR, and count value and functional analysis.
B_tbl <- read.table('Huebner201609_XENLAtx_table_B_clean.txt',sep = '\t', head = TRUE,row.names = 1) 
C_tbl <- read.table('Huebner201609_XENLAtx_table_C_clean.txt',sep = '\t', head = TRUE,row.names = 1)

GO_ready <- intersect(rownames(B_tbl),rownames(C_tbl))
write.table(GO_ready,'Huebner201609_XENLAtx_reference_gene.txt',quote=FALSE,row.names=FALSE,col.names = FALSE)

B_DMZ_sig_pos <- B_tbl[2:11][c(total_gene_pos),]
B_DMZ_sig_neg <- B_tbl[2:11][c(total_gene_neg),]
C_DMZ_sig_pos <- C_tbl[2:11][c(total_gene_pos),]
C_DMZ_sig_neg <- C_tbl[2:11][c(total_gene_neg),]
pos_value.table <- cbind(B_DMZ_sig_pos, B_sig_pos, C_DMZ_sig_pos, C_sig_pos)
neg_value.table <- cbind(B_DMZ_sig_neg, B_sig_neg, C_DMZ_sig_neg, C_sig_neg)
pos_value.table <- cbind(gene.names = row.names(pos_value.table), pos_value.table)
neg_value.table <- cbind(gene.names = row.names(neg_value.table), neg_value.table)
write.table(pos_value.table,"Huebner201609_XENLAtx_positive_final_table.txt",sep='\t',quote = FALSE,row.names = FALSE)
write.table(neg_value.table,"Huebner201609_XENLAtx_negative_final_table.txt",sep='\t',quote = FALSE,row.names = FALSE)
