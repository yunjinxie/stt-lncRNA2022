##=========================对hyper_test结果p矫正
setwd("J:/xieyunjin/test-GBM-model")
library(data.table)
ht<-fread("lnc-lnc-hyper-test.txt",sep="\t",data.table=F)
############### fdr校正 ###############
fdr_value<-p.adjust(ht[,6],"fdr")
########################################
hf<-cbind(ht,fdr_value)
########### fdr < 0.05 #################
hf_05<-hf[hf[,7]<0.05,]
num05<-nrow(hf_05)
lnc05<-c(as.character(hf_05[,1]),as.character(hf_05[,3]))
nlnc05<-length(unique(lnc05))
re1<-c("GBM",num05,nlnc05)#cancer lncRNA-lncRNA对数和lncRNA数
write.table(hf_05,"lnc-lnc-hyper-test05.txt",quote=F,sep="\t",col.names=F,row.names=F)

########## fdr < 0.01 ##################
hf_01<-hf[hf[,7]<0.01,]
num01<-nrow(hf_01)
lnc01<-c(as.character(hf_01[,1]),as.character(hf_01[,3]))
nlnc01<-length(unique(lnc01))
res1<-c("GBM",num01,nlnc01)#cancer lncRNA-lncRNA对数和lncRNA数
write.table(hf_01,"lnc-lnc-hyper-test01.txt",quote=F,sep="\t",col.names=F,row.names=F)

