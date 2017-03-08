library(ggplot2)
f_candi = read.table('candidate_gene_list.txt')
candi.list = sort(t(f_candi))
f_human_dat = read.table('core_human_gene_exp_listi_0224.txt', row.names= 2, header = TRUE)
f_human_dat = f_human_dat[c(-1)]
f_mouse_dat = read.csv('core_human_gene_exp_in_mouse_list_0224.csv',row.names= 2,header = TRUE,sep='\t')
f_mouse_dat = f_mouse_dat[c(-1)]
gene_list = rownames(f_human_dat)
sample_type = c()
sample_region = c()
sample_FC = c()
dat = c()
pval_bar = list()
asterisk = list()
for(i in 1:length(candi.list)){
  for(j in 1:length(gene_list)){
    if(candi.list[i] == gene_list[j]){
      sample.SN = log(c(as.numeric(f_human_dat[j,][1:6])),base = 2)
      sample.VTA = log(c(as.numeric(f_human_dat[j,][7:12])),base = 2)
      sample.val = max(c(sample.SN,sample.VTA))
      sample.SN.gene = paste(gene_list[j],'.SN',sep="")
      sample.VTA.gene = paste(gene_list[j],'.VTA',sep="")
      sample_type = c(sample_type,rep(c(sample.SN.gene,sample.VTA.gene),c(6,6)))
      sample_region = c(sample_region,rep(c('SN','VTA'),c(6,6)))
      if(median(sample.SN) > median(sample.VTA)){
        sample_FC = c(sample_FC,rep(c('Positive'),c(12)))
      }
      else if(median(sample.SN) < median(sample.VTA)){
        sample_FC = c(sample_FC,rep(c('Negative'),c(12)))
      }
      else{
        sample_FC = c(sample_FC,rep(c('NA'),c(12)))
      }
      dat = c(dat,sample.SN,sample.VTA)
      df <- data.frame(a=c(2*i-1,2*i-1,2*i,2*i),b= c(sample.val+0.2,sample.val+0.4,sample.val+0.4,sample.val+0.2))
      pval_bar[[i]] <- geom_line(data = df, aes(x=a, y=b))
      asterisk[[i]] <- annotate("text", x= 2*i-1+0.5, y= sample.val+0.5, label="*",size = 8)
      }
    }
}

temp = data.frame(val = dat,type=sample_type,region = sample_region, foldchange = sample_FC)
p <- ggplot(data = temp, aes(x=type, y=val)) +geom_boxplot(aes(fill=region),na.rm = TRUE)+ theme(text = element_text(size = 20),axis.text.x = element_text(angle=90, hjust = 1))+
 xlab("Gene")+ylab('log2 normalized value')+facet_grid(~foldchange,space="free",scales = "free")
plot(p)
ggsave('Human1.pdf',width = 20,height = 20)
