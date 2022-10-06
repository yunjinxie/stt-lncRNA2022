##======================寻找免疫协同大网中hub lncRNA的药靶基因=============================
##=========寻找靶基因=========
setwd("J:/xieyunjin/new_result/immune_hub/drug")
hubData<-read.table("hub_classify_syn_in_cancer.txt",sep='\t',header=T,as.is=T)
hub<-as.matrix(hubData$lncRNA)
cancer<-unique(unlist(strsplit(hubData[,2],",")))##在哪些癌症中协同
d<-"J:/xieyunjin/project_all/project_all/drug_target_result/lncRNA_gene"##所有功能对应的靶基因文件
f<-paste(cancer,"_lncRNA_gene.txt",sep='')
path<-file.path(d,f)
##药靶基因
drug<-read.table("J:/xieyunjin/project_all/project_all/drug_target_result/data/result_DrugBankID_HGNCsymbol_L01.txt",sep="\t",header=T,as.is=T)
colnames(drug)<-c("drug","target")
tarRes<-c()
for(i in 1:length(cancer)){
	geneData<-read.table(path[i],sep='\t',header=T,as.is=T)
	geneData<-as.matrix(geneData)
	allLnc<-substr(geneData[,"lncRNA"],1,15)
	hub2<-hub[grep(cancer[i],hubData[,2])]
	res<-lapply(hub2,function(x) cbind(x,unique(geneData[which(allLnc %in% x),"mRNA"])))
	res2<-do.call(rbind,res)##转化成矩阵
	tarRes<-rbind(tarRes,cbind(cancer[i],res2))
}
colnames(tarRes)<-c("cancer","lncRNA","target")
tarRes2<-merge(tarRes,hubData,by='lncRNA')
write.table(tarRes2,"immune_hub_target.txt",sep='\t',quote=F,row.names=F)
tarRes3<-unique(tarRes2[,c("lncRNA","target","type")])
drugRes<-merge(tarRes3,drug,by='target')
write.table(drugRes,"immune_hub_drugTarget.txt",sep='\t',quote=F,row.names=F)
######
##===========统计每一类hub调控靶基因个数=============
library(plyr)
count<-ddply(drugRes,.(type,lncRNA),function(x) length(unique(x$target)))
colnames(count)<-c("type","lncRNA","number")
write.table(count,"3type_hub_drugCount.txt",sep='\t',quote=F,row.names=F)
##===========画箱式图=====================
library(ggplot2)
library(ggpubr)
##顺序 common specific other!!!!!!!!
count<-read.table("3type_hub_drugCount.txt",sep='\t',header=T,as.is=T)
spe_count<-count[count[,1]=='specific',3]
com_count<-count[count[,1]=='common',3]
oth_count<-count[count[,1]=='other',3]
num<-list(spe_count,com_count,oth_count)
a<-matrix(c(1,2,2,3,1,3),ncol=2,byrow=T)
p<-apply(a,1,function(x) t.test(num[[x[1]]],num[[x[2]]])$p.value)
names(p)<-unlist(lapply(list(c("specific","common"),
                                   c("common","other"),
                                   c("specific","other")),function(x) paste(x[1],x[2],sep='_')))
write.table(p,"t_test.txt",sep='\t',quote=F,col.names=F)
####
col<-c("#607B8B","#FA8072","#8B7E66")
names(col)<-c("specific","common","other")
ggplot(count,aes(x=factor(type,levels=c("common","specific","other")),y=number,fill=type))+
      geom_boxplot(outlier.shape=NA)+
       theme_bw()+
      scale_fill_manual(values=col)+
      theme(legend.position="none")+
      labs(x=' ',y="number of drugs")+
      geom_signif(comparisons=list(c("specific","common"),
                                   c("common","other"),
                                   c("specific","other")), 
                  map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),test="t.test")+
      theme(plot.margin=unit(c(1,1,1,1),"cm"))
savePlot("immune_hub_drug","pdf",device=dev.cur(),restoreConsole=TRUE)
################
ggplot(count,aes(x=factor(type,levels=c("common","specific","other")),y=number,col=type))+
      geom_boxplot(outlier.shape=NA)+geom_jitter(shape=16, position=position_jitter(0.2))+
       theme_bw()+
      scale_colour_manual(values=col)+
      theme(legend.position="none")+
      labs(x=' ',y="number of drugs")+
      geom_signif(comparisons=list(c("specific","common"),
                                   c("common","other"),
                                   c("specific","other")), 
                  map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),test="t.test")+
      theme(plot.margin=unit(c(1,1,1,1),"cm"))
savePlot("immune_hub_drug2","pdf",device=dev.cur(),restoreConsole=TRUE)
################
pdf("immune_hub_drug_boxplot.pdf")
count$type<-factor(count$type,levels=c("common","specific","other"))
count[order(count$type),]
boxplot(number~type,count,col=c("common"="#FA8072","specific"="#607B8B","other"="#8B7E66"),
		las=1,ylab='number')
dev.off()
