#========================lnc2cancer3.0更新======================
setwd("J:/xieyunjin/new_result/immune_hub/lnc2cancer")
#3.0和文章搜的实验证实的
data <- readxl::read_excel("hub_cancer.xlsx")
data2 <- data[-1,c(1,3,6)]
data2 <- na.omit(data2)
colnames(data2) <- c("lncRNA","symbol","experi")

cla<-read.table("J:/xieyunjin/new_result/immune_hub/hallmark/data/hub_classify.txt",sep='\t',header=T,as.is=T)

common<-as.matrix(cla[cla$type=="common",])
specific<-as.matrix(cla[cla$type=="specific",])
other<-as.matrix(cla[cla$type=="other",])

data_cla <- merge(data2,cla,by="lncRNA")

com_dis<-as.matrix(data_cla[data_cla$type=="common",]$symbol)
spe_dis<-as.matrix(data_cla[data_cla$type=="specific",]$symbol)
oth_dis<-as.matrix(data_cla[data_cla$type=="other",]$symbol)

result<-c(com_dis,spe_dis,oth_dis)
type<-rep(c("common","specific","other"),c(length(com_dis),length(spe_dis),length(oth_dis)))
result2<-data.frame(result,type)
colnames(result2)<-c("lncRNA","type")

total<-c(nrow(common),nrow(specific),nrow(other))
p_data<-c(length(com_dis),length(spe_dis),length(oth_dis))/total
names(p_data)<-c("common","specific","other")

pdf("immune_hub_cancer_3.0.pdf")
barplot(p_data,col=c("common"="#FA8072","specific"="#607B8B","other"="#8B7E66"),las=1,ylab='percent')
text(0.25,seq(0.02,0.2,length.out = length(com_dis)),com_dis,cex=0.8,adj = 0)
text(1.45,seq(0.02,0.35,length.out = length(spe_dis)),spe_dis,cex=0.8,adj = 0)
text(2.7,seq(0.02,0.12,length.out = length(oth_dis)),oth_dis,cex=0.8,adj = 0)
dev.off()
