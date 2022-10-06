##==================免疫R1R2R3例子共享lncRNA表达==================
##共享lncRNA表达与其他lncRNA比较
library(data.table)
setwd("F:/xieyunjin/new_result/immune_R1R2R3/example/BLCA-CESC-HNSC-LUSC/BLCA-CESC-LUSC")
lnc<-read.table("lnc.txt",sep='\t',header=F,as.is=T)
lnc<-substr(as.matrix(lnc),1,15)
cancer<-c("BLCA","CESC","LUSC")
add<-"F:/xieyunjin/lncRNA_expression"
file<-list.files(add)
file2<-file[grep(paste(cancer,collapse="|"),file)]
fileAD<-file.path(add,file2)
add2<-"F:/xieyunjin/project_all/project_all/degree/cancer_degree/"
dir.create("pdf")

###
pdf(paste("pdf/other_expression.pdf",sep=''))
par(mfrow=c(1,3))
for(i in 1:length(cancer)){
  exp<-fread(fileAD[i],sep='\t',data.table=F)
  exp2<-as.matrix(exp[,-1])
  rownames(exp2)<-substr(as.matrix(exp[,1]),1,15)
  exp_lnc<-rownames(exp2)
  com_exp<-rowMeans(exp2[lnc,])
  other_lnc<-setdiff(exp_lnc,lnc)
  other_exp<-rowMeans(exp2[other_lnc,])
  p<-t.test(com_exp,other_exp,paired=F)$p.value
  boxplot(list(com_exp,other_exp),varwidth=T,las=1,col=c("#10C47E","gray"),ylab="the mean of expression(log2)",main=cancer[i],xaxt='n')
  axis(1,at=c(1,2),labels=c("common","other"))
}
dev.off()