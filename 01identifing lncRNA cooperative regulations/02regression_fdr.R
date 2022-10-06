###regression PÖµ¹ıÂË#####
setwd("J:/xieyunjin/test-GBM-model")
library(data.table)

#############
f<-paste("GBM_lm_1_lncRNA_result_part",1:3,".txt",sep='')
da_all<-c()
for(i in 1:length(f)){
	da<-fread(f[i],sep='\t',data.table=F)
	da_all<-rbind(da_all,da)
}
da_all2<-da_all[da_all[,5]<0.01,]
write.table(da_all2,"GBM_lncRNA_mRNA_fdr_2.txt",sep='\t',quote=F,row.names=F)

#############gene_id,lncRNA_id,lncRNA_coefficients,lncRNA_Pvalue,lncRNA_FDR
p<-0.1/dim(da_all2)[1]
da_all3<-da_all2[da_all2[,4]<p,1:3]
write.table(da_all3,"bon_2_GBM_lncRNA.txt",sep='\t',quote=F,row.names=F)