######### ��lncRNA��gene���й�ϵ��Ϊ�й�ϵ ##########################
######### ����\2013.9.2papers proms\task2\programs\regulator_target.R ###
######### ͬʱҲΪ�Ľ���
#ԭ�ļ�����һ��mRNA �ڶ���lncRNA ########################################
setwd("J:/xieyunjin/test-GBM-model")
library(data.table)
path1<-"J:/xieyunjin/project_all/project_all/bonforrni-lncRNA-0.1/bon_2_GBM_lncRNA.txt"

cf<-fread(path1,sep="\t",data.table=F)
f<-factor(cf[,2]) #lncRNA
subs<-split(cf[,1],f)#mRNA
NROW<-length(subs)
for(i in 1:NROW){
	r<-c(names(subs[i]),as.character(unique(unlist(subs[i]))))
	write.table(t(r),"GBM_convert_lnc_mRNA.txt",append=T,sep="\t",row.names=F,col.names=F,quote=F)
}
#�õ����ǵ�һ����lncRNA�������������صİ�gene