#===========================统计免疫泛癌网络中调控免疫功能的对数
library(plyr)
library(RColorBrewer)
display.brewer.all()

setwd("F:/xieyunjin/new_result/immune_type")

pan_go<-read.table("F:/xieyunjin/new_result/go_lnc_pairs_count/pan_cancer_lnc_go.txt",sep='\t',header=F,as.is=T)
colnames(pan_go)<-c("GO","go_names","pairs")
immu_go<-read.table("F:/xieyunjin/lncRNA_synergisty_project/immune_go/immune_go_term2.txt",sep='\t',header=T,as.is=T)

immu_pan_go<-merge(pan_go,immu_go,by="go_names")
data<-immu_pan_go[,c(1,3,5)]
write.table(data,"pan_immune_go.txt",sep='\t',quote=F,row.names=F)
#===========================统计泛癌网络中调控免疫功能的lncRNA个数
library(plyr)
library(RColorBrewer)
display.brewer.all()

setwd("F:/xieyunjin/new_result/immune_type")
data<-read.table("pan_immune_go.txt",sep='\t',header=T,as.is=T)
#######！！！！！！
result<-ddply(data,.(type,go_names),function(x) myCount(x))##计算调控每个免疫类型中免疫功能的lncRNA个数
colnames(result)<-c("type","go_names","number")
write.table(result,"immune_type_number_lnc.txt",sep='\t',quote=F,row.names=F)
##
myCount<-function(x){
  lnc<-apply(as.matrix(x$pairs),1,function(x) unlist(strsplit(x,"_")))
  len<-length(unique(c(lnc)))
  return(len)
}
#######画柱状图

#savePlot("immune_go","pdf",device=dev.cur(),restoreConsole=TRUE)
##box(col = 'black')
##tolower()大写改小写
##toupper()小写改大写
########14个免疫类型
result<-ddply(data,.(type),function(x) myCount(x))
colnames(result)<-c("type","number")
write.table(result,"type_number_lnc.txt",sep='\t',quote=F,row.names=F)
#############
ts<-read.table("type_number_lnc.txt",sep='\t',header=T,as.is=T)
#col<-colorRampPalette(c("green","blue","yellow","red","purple"))(5)
#col<-c(brewer.pal(8,"Set2")[1:7],brewer.pal(8,"Accent")[1:7])
col<-read.table("color.txt",sep='\t',header=F,as.is=T,comment.char="/")
colnames(col)<-c("color","type")
dataP2<-merge(ts,col,by='type')
label<-gsub("_"," ",dataP2[,1])
pdf("immune_type_number.pdf")
par(mai=c(1,4.5,1,.5))
barplot(dataP2[,2],horiz=T,col=dataP2[,3],las=1,names.arg=label,xlim=c(0,400))
dev.off()
#savePlot("immune_type","pdf",device=dev.cur(),restoreConsole=TRUE)
##
dataP<-read.table("immune_type_number_lnc.txt",sep='\t',header=T,as.is=T)
labels<-tolower(gsub("_"," ",dataP[,2]))
pdf("immune_go_names_number.pdf")
par(mai=c(1,5,1,.5))
barplot(dataP[,3],names.arg=labels,las=1,cex.names=.8,horiz=T,
        col=brewer.pal(9,"YlGn")[6],xlim=c(0,400))
dev.off()