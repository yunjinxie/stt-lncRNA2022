#===========任意两个癌症免疫网络中lncRNA做癌症间R1R2R3聚类=============
setwd("F:/xieyunjin/new_result/immune_R1R2R3/immune_lncRNA")
add<-"F:/xieyunjin/new_result/17_immune_network/node"
file<-dir(add)
cancer<-gsub("_node.txt","",file)
fileAD<-file.path(add,file)##癌症免疫网络节点
#======================================================================================================
##R1 任意两个癌症免疫网络中lncRNA交集并集
resultR1<-matrix(0,nrow=length(cancer),ncol=length(cancer))
rownames(resultR1)<-cancer
colnames(resultR1)<-cancer

resultR1min<-matrix(0,nrow=length(cancer),ncol=length(cancer))
rownames(resultR1min)<-cancer
colnames(resultR1min)<-cancer

resultR1max<-matrix(0,nrow=length(cancer),ncol=length(cancer))
rownames(resultR1max)<-cancer
colnames(resultR1max)<-cancer

for(i in 1:length(cancer)){
  nodeD<-read.table(fileAD[i],sep='\t',header=T,as.is=T)
  node<-as.matrix(nodeD[,1])
  for(j in 1:length(cancer)){
    nodeD2<-read.table(fileAD[j],sep='\t',header=T,as.is=T)
    node2<-as.matrix(nodeD2[,1])
    inter_node<-intersect(node,node2)
    union_node<-union(node,node2)
    resultR1[cancer[i],cancer[j]]<-length(inter_node)/length(union_node)
    resultR1max[cancer[i],cancer[j]]<-length(inter_node)/max(length(node),length(node2))
    resultR1min[cancer[i],cancer[j]]<-length(inter_node)/min(length(node),length(node2))
    
  }
}
write.table(resultR1,"R1_net.txt",sep='\t',quote=F)
write.table(resultR1max,"R1max_net.txt",sep='\t',quote=F)
write.table(resultR1min,"R1min_net.txt",sep='\t',quote=F)


#==========================================================================================================
immu_node<-read.table("immune_network_node.txt",sep='\t',header=F,as.is=T)
###R2 计算任意两个癌症的免疫lncRNA表达值相关性,用所有的免疫lncRNA
library(data.table)
result<-matrix(0,nrow=length(cancer),ncol=length(cancer))
rownames(result)<-cancer
colnames(result)<-cancer

###表达谱地址
exp_add<-"F:/xieyunjin/lncRNA_expression"
p<-list.files(exp_add)
p2<-p[grep(paste(cancer,collapse="|"),p)]
expAD<-file.path(exp_add,p2)

for(i in 1:length(cancer)){
  ##第一个癌症
  can<-fread(expAD[i],sep='\t',header=T,data.table=F)
  lnc<-substr(can[,1],1,15)#取lncRNA ID如 ENSG00000131484，去掉.之后的
  m1<-unlist(lapply(as.matrix(immu_node),function(x,lncRNA,exp) countMean(x,lnc,can)))
  ##第二个癌症，同上
  for(j in 1:length(cancer)){
    can2<-fread(expAD[j],sep='\t',header=T,data.table=F)  
    lnc2<-substr(can2[,1],1,15)
    m2<-unlist(lapply(as.matrix(immu_node),function(x,lncRNA,exp) countMean(x,lnc2,can2)))
    result[cancer[i],cancer[j]]<-cor.test(m1,m2)$estimate##计算任意基于大网的任意两个癌症的lncRNA表达值相关性
  }
}
write.table(result,"R2_net.txt",sep='\t',quote=F)

#####计算表达均值
countMean<-function(x,lncRNA,exp){##输入免疫lncRNA node和癌症表达谱中的lncRNA
  index<-grep(x,lncRNA,perl=T)
  if(length(index)!=0){
    m<-mean(as.numeric(as.matrix(exp[index,-1])))
    return(m)
  }else{
    return(0)
  }
}
#===============================================================================================================
###R3 计算R1 R2中值
f<-c("R1_net","R1min_net","R1max_net")
R2<-read.table("R2_net.txt",sep='\t',header=T,as.is=T)
for(e in f){
  f1<-read.table(paste(e,".txt",sep=''),sep='\t',header=T,as.is=T)
  R3<-(f1+R2)/2
  rownames(R3)<-rownames(f1)
  colnames(R3)<-colnames(f1)
  write.table(R3,paste(e,"_R3.txt",sep=''),sep='\t',quote=F)
  
}

##=========immune_R1R2R3画聚类图==========
dataP<-read.table("R1_net_R3.txt",sep='\t',header=T,as.is=T)
pdf("R1_net_R3.pdf")
pheatmap(dataP,col=c(colorRampPalette(c("#4575B4","#74ADD1","#FEE090"))(30),
                     colorRampPalette(c("#FEE090","#F46D43"))(30),colorRampPalette(c("#F46D43","#D73027"))(60)),
         clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",
         clustering_method="complete")