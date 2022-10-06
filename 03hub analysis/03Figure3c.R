#################统计hub lncRNA调控的GO及Hallmark及个数
setwd("J:/xieyunjin/new_result/immune_hub/hallmark")
hub<-read.table("data/immune_big_hub.txt",sep='\t',header=T)##lnc degree
hub<-as.matrix(hub)
hub_lnc<-hub[,1]
lnc_hall<-read.table("data/lnc_hallmark_goNames.txt",sep='\t',header=T,as.is=T)##lnc1 lnc2 hallmark go cancer
hall<-as.matrix(lnc_hall)
hall_lnc<-unique(c(hall[,1],hall[,2]))
inter_lnc<-intersect(hub_lnc,hall_lnc)
res1<-c()
for(i in 1:length(inter_lnc)){
    lnc<-inter_lnc[i]
    index<-which(hall[,1]==lnc|hall[,2]==lnc)
    go<-unique(hall[index,4])
    hallmark<-unique(hall[index,3])
    res<-cbind(lnc,paste(c(go),collapse=','),length(go),paste(c(hallmark),collapse=','),length(hallmark))
    res1<-rbind(res1,res)
}
colnames(res1)<-c("hub_lncRNA","go","go_len","hallmark","hallmark_len")
write.table(res1,"data/hub_immune_lnc_hall.txt",quote=F,row.names=F,sep='\t')
###########################统计3类hub分别调控的Hallmark个数
setwd("J:/xieyunjin/new_result/immune_hub/hallmark")
data<-read.table("data/hub_immune_lnc_hall.txt",sep='\t',header=T)
data<-as.matrix(data)
hall<-unique(unlist(strsplit(data[,4],",")))
hub<-read.table("data/hub_classify.txt",sep='\t',header=T)
common<-as.matrix(hub[hub$type=="common",1])
specific<-as.matrix(hub[hub$type=="specific",1])
other<-as.matrix(hub[hub$type=="other",1])

res<-matrix(0,nrow=length(hall),ncol=3)
for(i in 1:length(hall)){
  lnc<-data[grep(hall[i],data[,4]),1]
  common_hall_len<-length(intersect(common,lnc))
  specific_hall_len<-length(intersect(specific,lnc))
  other_hall_len<-length(intersect(other,lnc))
  res[i,1:3]<-c(common_hall_len,specific_hall_len,other_hall_len)
}

rownames(res)<-hall
colnames(res)<-c("common_hub","specific_hub","other_hub")

write.table(res,"data/immune_hallmark_hub.txt",sep='\t',quote=F)
###=============================计算3类hub占总共比例
per<-apply(res,1,function(x) x/sum(x))
colnames(per)<-gsub("_"," ",colnames(per))
dd<-c(length(common),length(specific),length(other))##总共识别的hub个数
ddP<-matrix(dd/sum(dd),ncol=1)
rownames(ddP)<-c("common_hub","specific_hub","other_hub")
colnames(ddP)<-"hub"
write.table(ddP,"data/all_hub_per.txt",sep='\t',quote=F)
##
per2<-cbind(per,ddP)
pdf("hallmark_hub.pdf")
par(mai=c(1,3,1,1))
barplot(per2,horiz=T,col=c("#FA8072","#607B8B","#8B7E66"),las=1)
dev.off()