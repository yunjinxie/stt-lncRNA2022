#=======================LUSC一致性聚类(km,euclidean)========================
setwd("J:/xieyunjin/new_result/ConsensusClusterPlus")
default <- getwd()
####对样本一致性聚类####
library(ConsensusClusterPlus)

SKCM_immune_lnc <- read.table("J:/xieyunjin/new_result/17_immune_network/node/SKCM_node.txt",sep='\t',header = T,as.is = T)
SKCM_immune_lnc <- as.matrix(SKCM_immune_lnc[,1])

exp <- read.table("J:/xieyunjin/lncRNA_expression/SKCM_regression_1_lncRNA.txt",sep = "\t",header = T,as.is = T)
exp2<-as.matrix(exp[,-1])
rownames(exp2) <- as.matrix(exp[,1])
exp_lnc<-rownames(exp2)

inter_lnc<-intersect(exp_lnc,SKCM_immune_lnc)
exp3<-exp2[inter_lnc,]
rs<-apply(exp3,1,function(x) length(which(x==-33)))##统计每个lncRNA样本值为0的个数
exp4<-exp3[rs<round(dim(exp3)[2]*0.5),]##阈值0.5

v<-apply(exp4,1,function(x) var(x))
v_sort<-sort(v,decreasing = T)
v_order<-order(v,decreasing = T)
cutoff<-round(length(v_order)*0.1)
while(v_sort[cutoff]==v_sort[cutoff+1]){
  cutoff<-cutoff+1
}
v_order2<-v_order[1:cutoff]
exp5<-exp4[v_order2,]

hub <- read.table("J:/xieyunjin/new_result/17_immune_network/hub2/SKCM_hub.txt",header = T,sep='\t',as.is = T)
hub <- as.matrix(hub[,1])
lnc <- intersect(hub,substr(rownames(exp5),1,15))
rownames(exp5) <- substr(rownames(exp5),1,15)
exp6 <- exp5[lnc,]

write.table(rownames(exp6),"SKCM_new/SKCM_filter_lnc.txt",sep='\t',quote = F,row.names = F,col.names = F)

title<-"SKCM_new"
results = ConsensusClusterPlus(exp6,maxK=5,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="euclidean",seed=1262118388.71279,plot="pdf")
setwd("SKCM_new")
for(j in 2:length(results)){
  write.table(results[[j]]$consensusClass,paste("cluster",j,".txt",sep=''),sep='\t',quote=F,col.names=F)
}
setwd(default)

#####生存####
library(survival)
library(survminer)
library(RColorBrewer)
clinic0<-read.table("J:/xieyunjin/UCSC_clinic_phenotype/TCGA-SKCM.survival.tsv",sep='\t',header = T,as.is = T)
clinic<-clinic0[,c(3,4,2)]
colnames(clinic)<-c("sample","time","status")
##两组
group<-read.table("SKCM_new/cluster2.txt",sep='\t',header=F,as.is=T)
sam<-apply(as.matrix(group[,1]),1,function(x) paste(unlist(strsplit(x,"\\."))[c(3,4,5)],collapse="-"))
group2<-data.frame(sample=sam,
                   group=group[,2])
clinic2<-merge(clinic,group2,by="sample")

table(clinic2$group)
#1   2 
#66 159 

km_fit<-survfit(Surv(time,status)~group,data=clinic2)
#km_p<-1-pchisq(survdiff(Surv(time,status)~group,data=clinic2)$chisq,1)

pdf("SKCM_new/SKCM_survival_table.pdf")
g<-ggsurvplot(km_fit,risk.table=T,#生存统计统计表
              # Change legends: title & labels
              legend.title = "Group",
              legend.labs = c("1","2"),
              conf.int=TRUE,#添加置信区间带
              palette = brewer.pal(9,"Set1")[c(3,5)],#颜色设置
              pval=TRUE,#log-rank检验
              pval.method=TRUE,#添加检验text
              pval.size=5,
              ncensor.plot=FALSE,# 画censor图
)+
  theme_survminer(font.legend=15,font.x=15,font.y=15,font.main=18)
print(g)
dev.off()
