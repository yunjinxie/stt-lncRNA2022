####========三类hub分别调控的免疫功能==========
library(plyr)
hc<-read.table("J:/xieyunjin/new_result/17_immune_network/bubble/data/hub_classify_syn_in_cancer.txt",sep='\t',header=T,as.is=T)
hc2<-hc[,c("lncRNA","type")]
hubGo<-read.table("hub_01_immuneType.txt",sep='\t',header=T,as.is=T)
hubGo2<-merge(hubGo,hc2)
write.table(hubGo2,"3hub_immune_type.txt",sep='\t',quote=F,row.names=F)

count<-ddply(hubGo2,.(immune_type,type),function(x) length(unique(x$lncRNA)))
colnames(count)<-c("immune_type","type","number")
write.table(count,"3hub_immune_type_number.txt",sep='\t',quote=F,row.names=F)
#####画玫瑰图
library(ggplot2)
count$type<-factor(count$type,levels=c("common","specific","other"))
count[order(count$type),]
ggplot(count,aes(x=immune_type,y=number,fill=type))+geom_bar(stat='identity')+coord_polar()+
  scale_fill_manual(values=c("specific"="#607B8B","common"="#FA8072","other"="#8B7E66"))+
  ylim(-10,80)+
  geom_hline(yintercept=0,col='gray')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0))
savePlot("rose","pdf",device=dev.cur(),restoreConsole=TRUE)


