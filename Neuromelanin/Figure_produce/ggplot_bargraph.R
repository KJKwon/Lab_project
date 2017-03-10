library(ggplot2)
f_candi = read.table('candidate_gene_list.txt')
candi.list = sort(t(f_candi))
f_human_dat = read.table('core_human_gene_exp_listi_0224.txt', row.names= 2, header = TRUE)
f_human_dat = f_human_dat[c(-1)]
f_mouse_dat = read.csv('core_human_gene_exp_in_mouse_list_0224.csv',row.names= 2,header = TRUE,sep='\t')
f_mouse_dat = f_mouse_dat[c(-1)]
gene.list = rownames(f_human_dat)
sample.list = c()
FC.list = c()
se = c()
for(i in 1:length(candi.list)){
  for(j in 1:length(gene.list)){
    if(candi.list[i] == gene.list[j]){
      sample.SN = c(as.numeric(f_human_dat[j,][13:15]))
      sample.VTA.md = median(c(as.numeric(f_human_dat[j,][16:18])))
      sample.SN_FC = sapply(sample.SN,function(x){log(x/sample.VTA.md,base=2)})
      FC.list = c(FC.list,mean(sample.SN_FC))
      sample.list = c(sample.list,gene.list[j])
      sample.sd.val = c(sd(sample.SN_FC)/sqrt(length(sample.SN_FC)))
      se = c(se,sample.sd.val)
    }
  }
}
temp = data.frame(gene = sample.list, val = FC.list,  se = se)
temp$colour = ifelse(temp$val < 0, "firebrick1","steelblue")
p = ggplot(temp,aes(x=gene,y=val, label = gene))+geom_bar(stat="identity",aes(fill=colour))+geom_errorbar(aes(ymin= val-se, ymax=val+se),width = 0.2,position = position_dodge(0.9))+
  coord_flip()
ggsave('Human2_bar.pdf',width = 20,height = 20)