##找到所有调控T细胞相关功能的lncRNA及所在癌症
setwd("F:/xieyunjin/new_result/B-Tcell/T-cell-diff-exp")

add<-"F:/xieyunjin/new_result/17_immune_network/network"
file<-dir(add)
file2<-file[grep("data",file)]
cancer<-gsub("_immune_data.txt","",file2)
fileAD<-file.path(add,file2)
cancer_names<-c()
for(i in 1:length(cancer)){
  can<-read.table(fileAD[i],sep='\t',header=T,as.is=T)
  lnc<-can[grep("T_CELL",can$go_names),c(2,3,5,7),drop=F]
  if(dim(lnc)[1]==0){
    next
  }
  output_dir<-paste("data",cancer[i],sep='/')
  dir.create(output_dir,recursive = T)
  setwd(output_dir)
  uni_lnc<-union(lnc[,1],lnc[,2])
  write.table(lnc,"T_network.txt",sep='\t',quote=F,row.names=F)
  write.table(uni_lnc,"lncRNA.txt",sep='\t',quote=F,row.names=F,col.names=F)
  cancer_names<-c(cancer_names,cancer[i])
  setwd("F:/xieyunjin/new_result/B-Tcell/T-cell-diff-exp")
  
}
write.table(cancer_names,"cancer.txt",sep='\t',quote=F,row.names=F,col.names=F)
#=======按照T(CD4,CD8)细胞比例做lncRNA差异表达分析，与调控T细胞激活功能的lncRNA做交集=============
####样本处理
library(data.table)
setwd("F:/xieyunjin/new_result/B-Tcell/T-cell-diff-exp")
cancer<-as.matrix(read.table("cancer.txt",sep='\t',header=F,as.is=T))
##表达谱
add<-"F:/xieyunjin/lncRNA_expression"
addAD<-paste(add,"/",cancer,"_regression_1_lncRNA.txt",sep='')
##免疫细胞比例
imEST<-read.table("F:/xieyunjin/new_result/immune_zhibiao_cor/data_processed/immuneEstimation.txt",sep='\t',header=T,as.is=T)
EST<-imEST[,c(2,3),drop=F]##样本  CD4 CD8
output_dir<-paste("data",cancer,sep='/')
sam<-rownames(EST)
##表达谱样本与细胞比例样本匹配

for(i in 1:length(cancer)){
  exp<-fread(addAD[i],sep='\t',data.table=F)
  sam1<-colnames(exp)[-1]
  sam2<-unlist(lapply(sam1,function(x) paste(unlist(strsplit(x,"\\."))[c(3,4,5)],collapse='-')))
  colnames(exp)<-c("gene_id",sam2)
  inter_sam<-intersect(sam,sam2)
  if(length(inter_sam)==0){
    next
  }
  exp2<-exp[,c("gene_id",inter_sam)]
  EST2<-EST[inter_sam,,drop=F]
  
  setwd(output_dir[i])
  write.table(exp2,"exp.txt",sep='\t',quote=F,row.names=F)
  write.table(EST2,"T_cell.txt",sep='\t',quote=F)
  setwd("F:/xieyunjin/new_result/B-Tcell/T-cell-diff-exp")
  
}

#=====将细胞免疫从小到大排序，按照中值分成两类，提取差异表达的lncRNA,
#============fold_change大于2且小于0.5，fdr<0.05
#============================与调控T细胞激活的lncRNA做交集===========================================
library(data.table)
library(pheatmap)
setwd("J:/xieyunjin/new_result/B-Tcell/T-cell-diff-exp")

cancer<-as.matrix(read.table("cancer.txt",sep='\t',header=F,as.is=T))
add<-paste("data",cancer,sep='/')
RPKM_add <- paste("J:/xieyunjin/RPKM_filter/",cancer,"_lncRNA_RPKM.txt",sep='')

output_dir<-paste("CD4/result",cancer,sep='/')
output_pdf<-"CD4/PDF"
dir.create(output_pdf)
count<-c()
for(i in 1:length(cancer)){
  BC<-read.table(paste(add[i],"/T_cell.txt",sep=''),sep='\t',header=T,as.is=T)
  BC2<-as.numeric(as.matrix(BC[,1]))#cd4
  cutoff<-median(BC2)
  if(cutoff>min(BC2)){
    group1<-rownames(BC)[BC2>=cutoff]
    group2<-rownames(BC)[BC2<cutoff]
  }else{
    group1<-rownames(BC)[BC2>cutoff]
    group2<-rownames(BC)[BC2<=cutoff]
  }
  #fold-change
  RPKM<-fread(RPKM_add[i],sep='\t',data.table=F,check.names = F)
  RPKM2<-as.matrix(RPKM[,-1])
  rownames(RPKM2)<-as.matrix(RPKM[,1])
  sam <- lapply(colnames(RPKM2),function(x) paste(unlist(strsplit(x,"\\."))[c(3:5)],collapse = "-"))
  colnames(RPKM2) <- sam

  Fold_change<-rowMeans(RPKM2[,group1])/rowMeans(RPKM2[,group2])
  #t.test
  exp <- fread(paste(add[i],"/exp.txt",sep=''),sep='\t',data.table=F)
  exp2<-as.matrix(exp[,-1])
  rownames(exp2)<-as.matrix(exp[,1])
  
  identical(rownames(RPKM2),rownames(exp2))
  #[1] TRUE
  
  t_pvalue<-apply(exp2,1,function(x) t.test(x[group1],x[group2],alternative = "two.side",conf.level = 0.95)$p.value)
  #fdr adjust
  fdr<-p.adjust(t_pvalue,method='fdr')
  
  result<-cbind("lnc"=names(t_pvalue),"p"=unname(t_pvalue),"fdr"=unname(fdr),"fold-change"=unname(Fold_change))
  dir.create(output_dir[i],recursive = T)
  setwd(output_dir[i])
  write.table(result,"result.txt",sep='\t',quote=F,row.names = F)
  #significant
  diff<-result[which(as.numeric(result[,4])>2|as.numeric(result[,4])<0.5),,drop=F]#fold-change
  diff_fdr<-diff[which(as.numeric(diff[,3])<0.05),,drop=F]#fdr BH
  write.table(diff_fdr,paste("exp_diff_result_fdr0.05.txt",sep=''),sep='\t',quote=F,row.names = F)
  
  setwd("J:/xieyunjin/new_result/B-Tcell/T-cell-diff-exp")
  
  #####=========================取交集 
  lnc<-read.table(paste(add[i],"lncRNA.txt",sep='/'),sep='\t',header=F,as.is=T)
  lnc <- as.matrix(lnc)
  
  inter<-intersect(lnc,diff_fdr[,1])
  syn_diff<-diff_fdr[which(diff_fdr[,1] %in% inter),]
  
  setwd(output_dir[i])
  write.table(syn_diff,"syn_CD4cell_diff.txt",sep='\t',quote=F,row.names=F)
  count<-rbind(count,cbind(cancer[i],nrow(diff_fdr),length(lnc),length(inter)))
  setwd("J:/xieyunjin/new_result/B-Tcell/T-cell-diff-exp")
  ####画聚类图
  if(length(inter)==0){
    next
  }
  dataP<-exp2[inter,c(group1,group2),drop=F]
  annotation_col=data.frame(
    type=factor(rep(c("group1","group2"),c(length(group1),length(group2))))
  )
  rownames(annotation_col)<-c(group1,group2)
  
  pdf(paste(output_pdf,"/",cancer[i],".pdf",sep=''))
  if(length(inter)==1){
    pheatmap(dataP,cluster_row=F,cluster_cols = F,
             col=c(colorRampPalette(c("blue","white"))(160),colorRampPalette(c("white","orange","red"))(50)),
             annotation_col = annotation_col)
  }else{
    pheatmap(dataP,cluster_cols = F,
             col=c(colorRampPalette(c("blue","white"))(160),colorRampPalette(c("white","orange","red"))(50)),
             annotation_col = annotation_col)
  }
  dev.off()
}
colnames(count)<-c("cancer","total_diff","Tcell_lnc","Tcell_lnc_diff")
write.table(count,"count_CD4cell.txt",sep='\t',quote=F,row.names=F)

##total_diff  Bcell_lnc Bcell_lnc_diff fisher精确检验，背景基因用全部基因组lncRNA基因，13870
library(ggplot2)
setwd("J:/xieyunjin/new_result/B-Tcell/T-cell-diff-exp")
cell<-c("CD4","CD8")
filenames<-paste("count_",cell,"cell",sep='')
dir.create("ORpdf")

result <- c()
for(i in 1:length(cell)){
  count<-read.table(paste(filenames[i],".txt",sep=''),sep='\t',header=T,as.is=T)
  count2<-as.matrix(count[,-1])
  rownames(count2)<-count[,1]
  p<-apply(count2,1,function(x) fisherP(x))
  OR<-apply(count2,1,function(x) fisherOR(x))
  lower<-apply(count2,1,function(x) fisherConf.lower(x))
  upper<-apply(count2,1,function(x) fisherConf.upper(x))
  result<-rbind(result,data.frame("cancer" = names(p),"cell" = cell[i],p,OR,lower,upper))
}

write.table(result,"T_fisher-exact.txt",sep='\t',quote=F,row.names = F)
#plot
result <- read.table("T_fisher-exact.txt",sep='\t',header=T,as.is = T)

result_sig<-result[result$p<0.05,]
result_sig$asterisk <- "*"
result_sig[result_sig$p <= 0.01&result_sig$p > 0.001,]$asterisk <- "**"
result_sig[result_sig$p < 0.001,]$asterisk <- "***"

#level <- sort(levels(result_sig$cancer),decreasing = T)
level <- unique(result_sig$cancer)

pdf("T_OR-2.pdf",width = 12,height = 5)
ggplot(result_sig,aes(x=log2(OR),y=factor(cancer,levels = level),colour = cell))+
  geom_point(cex=2)+
  geom_errorbarh(aes(xmin=log2(lower),xmax=log2(upper),height=.2,colour = cell))+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size=19,angle=90,colour = "black"),
        axis.text.y = element_text(size=15,colour = "black"),
        axis.title.y = element_text(size=19,colour = "black"),
        legend.text = element_text(size=14,colour="black"),
        legend.title = element_blank())+
  scale_colour_manual(values = c("CD4" = "#88c0c6","CD8" = "#bce1e5"))+
  labs(x="log2(OR)",y="",colour = "Cell")+
  geom_text(aes(label=asterisk,colour = cell))


dev.off()
#######function fisher.test
myfisher<-function(df){
  N<-13870
  a<-df[3]
  c<-df[2]-df[3]
  b<-df[1]-df[3]
  cancer<-names(df)
  td<-data.frame(gene.in.interest=c(a,c),gene.not.interest=c(b,N-a-b-c))
  res<-fisher.test(td)
  return(res)
}
#
fisherP<-function(df){
  fisher_result<-myfisher(df)
  p<-fisher_result$p.value
  return(p)
}
fisherOR<-function(df){
  fisher_result<-myfisher(df)
  OR<-fisher_result$estimate
  return(OR)
}
fisherConf.lower<-function(df){
  fisher_result<-myfisher(df)
  lower<-fisher_result$conf.int[1]
  return(lower)
}
fisherConf.upper<-function(df){
  fisher_result<-myfisher(df)
  upper<-fisher_result$conf.int[2]
  return(upper)
}


