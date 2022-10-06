####计算免疫指标的差异####
library(pheatmap)
library(RColorBrewer)

setwd("J:/xieyunjin/new_result/ConsensusClusterPlus/SKCM_new")

##免疫signatures
can_sam<-read.table("cluster2.txt",sep='\t',header=F,as.is=T)
sam<-apply(as.matrix(can_sam[,1]),1,function(x) paste(unlist(strsplit(x,"\\."))[c(3,4,5)],collapse='-'))
can_sam2<-cbind(sam,can_sam[,2])
colnames(can_sam2)<-c("sample","group")

measure0<-read.table("zhibiao_hebing.txt",sep='\t',header=T,as.is=T)
measure0 <- measure0[,-c(1,2)]
measure0_1 <- cbind(rownames(measure0),measure0)
colnames(measure0_1)[1] <- "sample"

check_point <- c("PDCD1","CD274","CTLA4")
m_exp <- read.table("J:/xieyunjin/regression_1_mRNA_data/expression/SKCM_regression_1_mRNA.txt",sep='\t',header = T,as.is = T,check.names = F)
m_exp[1:3,1:3]
m_exp2 <- m_exp[,-1]
rownames(m_exp2) <- as.matrix(m_exp[,1])
m_exp_log <- log2(m_exp2+0.001)
m_exp_c <- m_exp_log[check_point,]
colnames(m_exp_c) <- gsub("\\.","-",colnames(m_exp_c))
m_exp_c2 <- t(m_exp_c)
m_exp_c3 <- cbind.data.frame(rownames(m_exp_c2),m_exp_c2)
colnames(m_exp_c3)[1] <- "sample"

measure0_all <- merge(measure0_1,m_exp_c3,by= "sample")

measure1<-merge(measure0_all,can_sam2,by='sample')
measure1$group <- factor(measure1$group,levels = c(1,2))
dataP <- measure1[order(measure1$group),]

write.table(dataP,"SKCM_zhibiao_hebing2_plot.txt",sep='\t',quote=F,row.names=F)
result<-apply(dataP[,2:8],2,function(x) wilcox.test(x~dataP$group)$p.value)
result2<-data.frame("zhibiao"=names(result),"pvalue"=unname(result))
write.table(result2,"SKCM_zhibiao_hebing2_wilcox.txt",sep='\t',quote=F,row.names=F)

######
dataP2<-dataP[,c(6:8,3:5)]
#dataP2_log <- log2(dataP2)
rownames(dataP2)<-dataP[,1]
annotation_row=data.frame(Group=factor(as.character(dataP$group)))
rownames(annotation_row)<-dataP$sample
ann_colors=list(Group=c('1'=brewer.pal(9,"Set1")[3],'2'=brewer.pal(9,"Set1")[5]))
##
pdf("SKCM_zhibiao_hebing.pdf")
s<-c(seq(-1420,-11,by=10),seq(-10,-0.1,by=0.01),seq(0,2.5,by=0.01),seq(2.6,22,by=0.01),seq(23,150,by=0.01),seq(151,3270,by=10))
# col2<-c(colorRampPalette(c("#0000ff","#673dff"))(length(seq(-1420,-4,by=10))),
#         colorRampPalette(c("#673dff","#9265ff","#b38bff","#FFFFFF"))(length(seq(-3.9,-0.1,by=0.01))),
#         colorRampPalette(c("#FFFFFF","#ffdfd4","#ffbfaa"))(length(seq(0,2.5,by=0.01))),
#         colorRampPalette(c("#ffbfaa","#ff9e81","#ff7b5a"))(length(seq(2.6,22,by=0.01))),
#         colorRampPalette(c("#ff7b5a","#ff5232","#ff0000"))(length(seq(23,150,by=0.01))),
#         colorRampPalette(c("#ff0000","darkred"))(length(seq(151,3270,by=10))))

col2<-c(colorRampPalette(c("#1e90ff","#64a1ff","#8cb3ff"))(length(seq(-1420,-11,by=10))),
        colorRampPalette(c("#8cb3ff","#c9d8ff","#FFFFFF"))(length(seq(-10,0,by=0.01))),
        colorRampPalette(c("#FFFFFF","#fffacd"))(length(seq(0,2.5,by=0.01))),
        colorRampPalette(c("#fffacd","#ffbda5","#d67365"))(length(seq(2.6,22,by=0.01))),
        colorRampPalette(c("#d67365","#c54f43"))(length(seq(23,150,by=0.01))),
        colorRampPalette(c("#c54f43","#b22222"))(length(seq(151,3270,by=10))))

pheatmap(dataP2,annotation_row=annotation_row,cluster_cols = F,cluster_rows = F,show_rownames = F,breaks = s,
         col=col2,annotation_colors =  ann_colors)
dev.off()
################考虑肿瘤突变负荷###########
##免疫signatures
can_sam<-read.table("cluster2.txt",sep='\t',header=F,as.is=T)
sam<-apply(as.matrix(can_sam[,1]),1,function(x) paste(unlist(strsplit(x,"\\."))[c(3,4,5)],collapse='-'))
can_sam2<-cbind(sam,can_sam[,2])
colnames(can_sam2)<-c("sample","group")

measure0<-read.table("zhibiao_hebing.txt",sep='\t',header=T,as.is=T)
#measure0 <- measure0[,-c(1,2)]
measure0_1 <- cbind(rownames(measure0),measure0)
colnames(measure0_1)[1] <- "sample"

check_point <- c("PDCD1","CD274","CTLA4")
m_exp <- read.table("J:/xieyunjin/regression_1_mRNA_data/expression/SKCM_regression_1_mRNA.txt",sep='\t',header = T,as.is = T,check.names = F)
m_exp[1:3,1:3]
m_exp2 <- m_exp[,-1]
rownames(m_exp2) <- as.matrix(m_exp[,1])
m_exp_log <- log2(m_exp2+0.001)
m_exp_c <- m_exp_log[check_point,]
colnames(m_exp_c) <- gsub("\\.","-",colnames(m_exp_c))
m_exp_c2 <- t(m_exp_c)
m_exp_c3 <- cbind.data.frame(rownames(m_exp_c2),m_exp_c2)
colnames(m_exp_c3)[1] <- "sample"

measure0_all <- merge(measure0_1,m_exp_c3,by= "sample")

measure1<-merge(measure0_all,can_sam2,by='sample')
measure1$group <- factor(measure1$group,levels = c(1,2))
dataP <- measure1[order(measure1$group),]

write.table(dataP,"SKCM_zhibiao_hebing2_plot.txt",sep='\t',quote=F,row.names=F)
result<-apply(dataP[,2:10],2,function(x) wilcox.test(x~dataP$group)$p.value)
result2<-data.frame("zhibiao"=names(result),"pvalue"=unname(result))
write.table(result2,"SKCM_zhibiao_hebing2_wilcox_all.txt",sep='\t',quote=F,row.names=F)
