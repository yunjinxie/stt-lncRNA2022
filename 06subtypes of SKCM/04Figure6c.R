####计算两类样本间免疫细胞比例差异####
library(pheatmap)
library(ggplot2)
library(car)
library(RColorBrewer)
library(ggpubr)

setwd("J:/xieyunjin/new_result/ConsensusClusterPlus")

cell_type<-read.table("J:/xieyunjin/new_result/immune_zhibiao_cor/data_processed/immuneEstimation.txt",sep='\t',header=T,as.is=T)
cell_sam<-rownames(cell_type)

can_sam<-read.table("SKCM_new/cluster2.txt",sep='\t',header=F,as.is=T)
can_sam2<-as.matrix(can_sam[,-1])
sam<-apply(as.matrix(can_sam[,1]),1,function(x) paste(unlist(strsplit(x,"\\."))[c(3,4,5)],collapse='-'))
rownames(can_sam2)<-sam

cell_inter<-intersect(sam,cell_sam)
cell_type2<-cell_type[cell_inter,]
can_sam3<-can_sam2[cell_inter,,drop=F]

data<-data.frame("label"=factor(can_sam3),cell_type2)
write.table(data,"SKCM_new/SKCM_immune_cell_porportion.txt",sep='\t',quote = F)

dataP<-data.frame("label"=rep(data$label,(dim(data)[2]-1)),
                  "cell_type"=rep(colnames(data)[-1],each=dim(data)[1]),
                  "number"=as.numeric(as.matrix(data[,-1])))
dataP$group<-paste(dataP$label,dataP$cell_type,sep='_')
len<-dim(data)[2]-1
a<-paste(rep(c(1,2),times=len),rep(colnames(data)[-1],each=2),sep='_')
col_res<-c()
for(j in 1:length(a)){
  if(grepl("^1_",a[j])){
    col_res<-rbind(col_res,cbind(a[j],brewer.pal(9,"Set1")[3]))
  }else{col_res<-rbind(col_res,cbind(a[j],brewer.pal(9,"Set1")[5]))}
}
col_res2<-col_res[,2]
names(col_res2)<-col_res[,1]
com<-matrix(a,ncol=2,byrow=T)
p<-c()
for(k in 1:dim(com)[1]){
  n<-dataP[dataP$group==com[k,1],"number"]
  m<-dataP[dataP$group==com[k,2],"number"]
  
  p<-c(p,wilcox.test(n,m)$p.value)
}
res<-cbind("SKCM",com,p)
pdf("SKCM_new/SKCM_cell_type.pdf")

p1<-ggplot(dataP,aes(x=factor(group,levels=col_res[,1]),y=number,col=group))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4)+
  scale_colour_manual(values =col_res2)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60,vjust = 0.5),
        panel.grid = element_blank())+
  geom_signif(comparisons=list(com[1,],com[2,],com[3,],com[4,],com[5,],com[6,]),
              test="wilcox.test",map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05),y_position=rep(1,6))+
  ylim(0,1)+labs(x=" ",y="Proportion of cell")+
  theme(legend.position='none')+
  scale_x_discrete(labels=c("B cell","","CD4 T cell","","CD8 T cell","","Neutrophil",""," Macrophage","","Dendritic",""))+
  coord_fixed (15)

print(p1)
dev.off()

write.table(res,"SKCM_new/wilcox_result.txt",sep='\t',quote=F,row.names=F,col.names=F)
