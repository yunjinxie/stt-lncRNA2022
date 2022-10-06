#######################################
#超几何检验
setwd("J:/xieyunjin/test-GBM-model")
library(data.table)
path1<-"J:/xieyunjin/project_all/project_all/bonforrni-lncRNA-0.1/bon_2_GBM_lncRNA.txt"
cf<-fread(path1,sep='\t',data.table=F)
##癌症中mRNA总个数
mRNA<-unique(cf[,1])
N<-length(mRNA)
##===============================================
linc_gene<-readLines("GBM_convert_lnc_mRNA.txt")
lg<-strsplit(linc_gene,"\t")
lglen<-length(lg)#lncRNA numbers
LEN<-lglen-1
for(i in 1:LEN){
	for(j in (i+1):lglen){
		################# 超几何参数 #################
		m<-length(lg[[i]])-1
		n<-N-m+1
		k<-length(lg[[j]])-1
		x<-length(intersect(lg[[i]][-1],lg[[j]][-1]))
		################ 求p值 ######################
		if(x>1){
			p<-1-phyper(x,m,n,k)
			########### 结果(linc_i,i基因集基因数，linc_j,j基因集基因数，交集基因数，检验p_value) ###
			r<-c(lg[[i]][1],m,lg[[j]][1],k,x,p)
			write.table(t(r),"lnc-lnc-hyper-test.txt",append=T,col.names=F,row.names=F,quote=F,sep="\t")
		}
	}
}
