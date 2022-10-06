#===============================================17个癌症调控免疫类型lncRNA对数统计=============================================
##================行为调控免疫功能的17个癌症，列为14个免疫类型，矩阵代表在某个癌症中调控某个免疫类型的lncRNA对数=================
library(plyr)
setwd("F:/xieyunjin/new_result/17_immune_network")
immu_data<-read.table("F:/xieyunjin/lncRNA_synergisty_project/immune_go/immune_go_term2.txt",sep='\t',header=T,as.is=T)
immu_type<-unique(as.matrix(immu_data$type))

net_data<-dir("network")
net<-net_data[grep("data",net_data)]
cancer<-gsub("_immune_data.txt","",net)
netAD<-file.path("network",net)

result_pairs<-matrix(0,nrow=length(immu_type),ncol=length(net))
rownames(result_pairs)<-immu_type
colnames(result_pairs)<-cancer
title<-cbind("type","cancer","pairs")
write.table(title,"cancer_type_pairs.txt",quote=F,sep='\t',row.names=F,append=T,col.names=F)

for(i in 1:length(cancer)){
  can<-read.table(netAD[i],sep='\t',header=T,as.is=T)
  res<-ddply(can,.(type),function(x) data.frame(cancer[i],unique(x$pairs)))
  colnames(res)<-c("type","cancer","pairs")
  write.table(res,"cancer_type_pairs.txt",quote=F,sep='\t',row.names=F,append=T,col.names=F)
  count<-ddply(res,.(cancer,type),nrow)
  result_pairs[count[,2],i]<-count[,3]
}
write.table(result_pairs,"count_result_pairs_new.txt",quote=F,sep='\t')
#cs<-colSums(result)
##=========================画色图加柱状图=================================
library(pheatmap)
library(RColorBrewer)
display.brewer.all()
#per_result
count<-read.table("count_result_pairs_new.txt",header=T,sep='\t',as.is=T)
rs<-rowSums(count)
cs<-colSums(count)
count2<-count[order(rs,decreasing = T),order(cs,decreasing = T)]

sNum<-read.table("F:/xieyunjin/lncRNA_synergisty_project/cancer_network_count.txt",sep='\t',as.is=T,header=T)
sNum2<-sNum[-c(9,14,16),2]##对数
per_result<-t(t(count2)/sNum2[order(cs,decreasing = T)])
rownames(per_result)<-gsub("_"," ",rownames(per_result))
write.table(per_result,"per_count_result_sum(network)_order_pairs.txt",quote=F,sep='\t')

pdf("count_per_sum(network)_order.pdf",width=10,height=10)
pheatmap(per_result,cellwidth=20,cellheight=20,display_numbers=count2,cluster_cols = F,cluster_rows = F,
         col=c(colorRampPalette(c(brewer.pal(9,"Blues")[c(2,6)]))(50),colorRampPalette(c(brewer.pal(9,"Blues")[6],"orange","red"))(300)))
dev.off()
#======================分母是免疫协同网络中关系对====================
library(plyr)
setwd("G:/lncRNA协同调控/new_result/17_immune_network")

library(pheatmap)
library(RColorBrewer)
display.brewer.all()
#per_result
count<-read.table("juzhen/count_result_pairs_new.txt",header=T,sep='\t',as.is=T)
rs<-rowSums(count)
cs<-colSums(count)
count2<-count[order(rs,decreasing = T),order(cs,decreasing = T)]
cancer<-colnames(count2)

sNum<-read.table("count/count.txt",sep='\t',as.is=T,header=T)

per_result<-t(t(count2)/sNum[cancer,1])
rownames(per_result)<-gsub("_"," ",rownames(per_result))
write.table(per_result,"juzhen/per_count_result_sum(immune_network)_order_pairs.txt",quote=F,sep='\t')

pdf("juzhen/count_per_sum(immune_network)_order.pdf",width=10,height=10)
pheatmap(per_result,cellwidth=20,cellheight=20,fontsize = 13,display_numbers=count2,cluster_cols = F,cluster_rows = F,
         col=c(colorRampPalette(c(brewer.pal(9,"Blues")[c(2,6)]))(50),colorRampPalette(c(brewer.pal(9,"Blues")[6],"orange","red"))(300)))
dev.off()

###行：调控该功能的所有lncRNA对数
data<-read.table("cancer_type_pairs.txt",sep='\t',as.is=T,header=T)
rs<-ddply(data,.(type),function(x) length(unique(x$pairs)))
rs2<-rs[order(rs[,2],decreasing = F),]
pdf("type.pdf")
par(mai=c(1,4.5,1,0.1))
barplot(t(rs2[,2]),horiz=T,names.arg=rs2[,1],las=1,col=brewer.pal(9,"Blues")[4],border='gray',xlim=c(0,1000))
dev.off()
##列：在一个癌症中所有功能被调控几次
data<-read.table("juzhen/cancer_type_pairs.txt",sep='\t',as.is=T,header=T)
cs<-ddply(data,.(cancer),function(x) length(unique(x$pairs)))
cs2<-cs[order(cs[,2],decreasing = T),]
pdf("juzhen/cancer.pdf")
barplot(t(cs2[,2]),las=2,col=brewer.pal(9,"Blues")[4],border='gray')
dev.off()
