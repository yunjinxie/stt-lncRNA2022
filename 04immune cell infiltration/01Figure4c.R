#=============计算协同lncRNA表达和免疫细胞浸润比例相关性=============================
library(data.table)
library(psych)

setwd("J:/xieyunjin/new_result/immune_zhibiao_cor")
mea<-read.table("J:/xieyunjin/new_result/immune_zhibiao_cor/data_processed/immuneEstimation.txt",sep='\t',header=T,as.is=T)##免疫细胞类型
#表达谱
exp_data<-"J:/xieyunjin/lncRNA_expression"
exp_file<-dir(exp_data)
expAD<-file.path(exp_data,exp_file)
cancer<-gsub("_regression_1_lncRNA.txt","",exp_file)
#癌症协同网络
d<-"J:/xieyunjin/project_all/project_all/lncRNA_FDR_data/cancer_fit_model_p_adjust"
file<-dir(d)
path<-file.path(d,file)

output_dir<-"cooperative_cell_type"
dir.create(output_dir)
default<-getwd()

for(j in 1:length(cancer)){
  exp<-fread(expAD[j],sep='\t',data.table=F)
  exp2<-as.matrix(exp[,-1])
  rownames(exp2)<-as.matrix(exp[,1])
  
  net<-fread(path[j],sep=',',data.table = F)
  net<-as.matrix(net)
  lnc<-union(net[,1],net[,2])
  
  sam2<-as.matrix(colnames(exp2))
  sam22<-unlist(lapply(sam2,function(x) paste(unlist(strsplit(x,"\\."))[-c(1,2)],collapse='-')))
  colnames(exp2)<-sam22
  sam<-rownames(mea)
  inter_sam<-intersect(sam,sam22)
  if(length(inter_sam)==0){
    next
  }
  mea2<-mea[inter_sam,]
  exp3<-exp2[lnc,inter_sam]
  
  Pvalue<-matrix(0,nrow=length(lnc),ncol=dim(mea2)[2])
  estimate2<-matrix(0,nrow=length(lnc),ncol=dim(mea2)[2])
  ##默认的校正方法holm
  res<-corr.test(t(exp3), mea2, method= "spearman",adjust="holm")
  estimate<-res$r
  pvalue<-res$p
  
  setwd(output_dir)
  write.table(pvalue,paste(cancer[j],"_pvalue.txt",sep=''),sep='\t',quote=F)
  write.table(estimate,paste(cancer[j],"_estimate.txt",sep=''),sep='\t',quote=F)
  #===============阈值|R| > 0.3 and P < 0.05的显示相关系数，否则为0==============
  pvalue2<-as.numeric(pvalue)
  estimate2<-as.numeric(estimate)
  estimate2[!(pvalue2<0.05 & abs(estimate2)>0.3) ]<-0
  result<-matrix(estimate2,nrow=nrow(estimate))
  colnames(result)<-colnames(mea2)
  rownames(result)<-lnc
  write.table(result,paste(cancer[j],"_pvalue_absR3_05.txt",sep=''),sep='\t',quote=F)
  setwd(default)
}

####统计每个癌症中与免疫细胞浸润比例显著相关的lncRNA个数
setwd("J:/xieyunjin/new_result/immune_zhibiao_cor/cooperative_cell_type")

d<-"J:/xieyunjin/new_result/immune_zhibiao_cor/cooperative_cell_type"
file<-dir(d)
file2<-file[grep("abs",file)]
path<-file.path(d,file2)

output<-"count_result"
dir.create(output)
default<-getwd()

cancer<-gsub("_pvalue_absR3_05.txt","",file2)

count<-c()
por<-c()
for(i in 1:length(file2)){
  es<-read.table(path[i],sep='\t',header=T,as.is=T)
  len<-nrow(es)
  cat(len)
  cat("\n")
  count0<-apply(es,2,function(x) length(which(x!=0)))
  por0<-count0/len
  por<-rbind(por,por0)
  count<-rbind(count,count0)
}
rownames(count)<-cancer
rownames(por)<-cancer

setwd(output)
write.table(count,"count.txt",sep='\t',quote=F)
write.table(por,"porpotion.txt",sep='\t',quote=F)
setwd(default)
#==============================画图==============================
setwd("J:/xieyunjin/new_result/immune_zhibiao_cor/cooperative_cell_type/count_result")

library(corrplot)
library(RColorBrewer)
display.brewer.all()

por <- read.table("porpotion.txt",sep='\t',header = T,as.is = T)
por <- as.matrix(por)


pdf("porpotion.pdf")

col2 <- colorRampPalette(c("#F4A582","#FDDBC7", "#FFFFFF", "#D1E5F0",
                           "#92C5DE","#4393C3", "#2166AC", "#053061"))

colnames(por)<-gsub("_"," ",colnames(por))
corrplot(por,is.corr=FALSE,method="circle",
         col=col2(30),bg=brewer.pal(9,"Pastel1")[9],outline=F,cl.pos='r',
         cl.lim = c(0,1),cl.length =5,cl.cex = 1,cl.ratio = 0.5 ,cl.align = "r",
         tl.cex=1,tl.col='black')

dev.off()


