##全部样本画生存曲线，并km检验
library(data.table)
library(survival)
library(RColorBrewer)
library(survminer)

setwd("J:/xieyunjin/new_result/immune_R1R2R3/example/BLCA-CESC-HNSC-LUSC/BLCA-CESC-LUSC")
lnc<-read.table("lnc.txt",sep='\t',header=F,as.is=T)
lnc<-as.matrix(lnc)
##表达
cancer<-c("BLCA","CESC","LUSC")
add<-"J:/xieyunjin/new_result/single_lnc_cox/expression/"
expAD<-paste(add,cancer,"_expression.txt",sep='')
##生存
add2<-"J:/xieyunjin/UCSC_clinic_phenotype/"
surAD<-paste0(add2,"TCGA-",cancer,".survival.tsv",sep='')

cox_result2_all <- c()
for(i in 1:length(cancer)){
  exp<-fread(expAD[i],sep='\t',data.table=F)
  exp2<-as.matrix(exp[,-1])
  rownames(exp2)<-exp[,1]
  sur0<-fread(surAD[i],sep='\t',data.table=F)
  sur<-unique(sur0[,c(3,4,2)])
  colnames(sur)<-c("sample","time","status")
  ff <- table(sur$sample)
  max(ff)
  inter_sam <- intersect(sur$sample,colnames(exp2))
  exp3 <- exp2[,inter_sam]
  sur2 <- sur[match(inter_sam,sur$sample),]
  cox_result<-apply(lnc,1,function(x) singleCox(x,exp2 = exp3,sur = sur2))
  cox_result2<-t(cox_result)
  cox_result2_all <- rbind(cox_result2_all,cbind(cancer[i],cox_result2))
}
colnames(cox_result2_all) <- c("cancer","lncRNA","p_value","coef","HR","lower","upper")
write.table(cox_result2_all,"cox_result_lncRNA.txt",sep='\t',quote = F,row.names = F)
##
singleCox<-function(x,exp2,sur){
  exp3<-exp2[x,]
  data<-list(time=sur$time,status=sur$status,exp=exp3)
  cox_res<-summary(coxph(Surv(time,status)~exp,data))
  p_value<-cox_res$coefficients[,5]
  Coef<-cox_res$coefficients[,1]
  HR <- cox_res$coefficients[,2]
  lower <- cox_res$conf.int[,3]
  upper <- cox_res$conf.int[,4]
  result<-cbind(x,p_value,Coef,HR,lower,upper)
  return(result)
}

cox_result2_all <- read.table("cox_result_lncRNA.txt",sep='\t',header = T,as.is = T)
cox_result2_all_sig <- cox_result2_all[as.numeric(cox_result2_all[,3]) < 0.05,]

#########按表达分组(确定最佳分组)############
dir.create("bestSeparation")
cancer <- unique(cox_result2_all_sig[,1])
km_p_all <- c()
for(i in 1:length(cancer)){
  lnc <- cox_result2_all_sig[cox_result2_all_sig[,1] == cancer[i],2]

  add<-"J:/xieyunjin/new_result/single_lnc_cox/expression/"
  expAD<-paste(add,cancer[i],"_expression.txt",sep='')
  exp<-fread(expAD,sep='\t',data.table=F)
  exp2<-as.matrix(exp[,-1])
  rownames(exp2)<-exp[,1]
  
  add2<-"J:/xieyunjin/UCSC_clinic_phenotype/"
  surAD<-paste0(add2,"TCGA-",cancer[i],".survival.tsv",sep='')
  sur0<-fread(surAD,sep='\t',data.table=F)
  sur<-unique(sur0[,c(3,4,2)])
  colnames(sur)<-c("sample","time","status")
  
  for(j in 1:length(lnc)){
   
    lnc_exp <- exp2[lnc[j],,drop=F]
    lnc_exp2 <- cbind(sample = colnames(lnc_exp),t(lnc_exp))
    
    svdata <- merge(lnc_exp2,sur,by= "sample")
    ##对数据集的基因进行bestSeparation统计
    svdata[,2] <-as.numeric(svdata[,2])
    
    res.cut <- surv_cutpoint(svdata, time = "time",
                             event = "status",
                             variables = colnames(svdata)[2],
                             minprop = 0.3) #默认组内sample不能低于30%
    res.cat <- surv_categorize(res.cut)
    colnames(res.cat)[3] <- "group"
    res.cat2 <- cbind(svdata,res.cat)
    
    write.table(res.cat2,paste0("bestSeparation/",cancer[i],"_",lnc[j],"_bestSeparation.txt"),
                sep='\t',quote = F,row.names = F)
    
    kmfit<-survfit(Surv(time,status)~group,res.cat)
    km_p<-1-pchisq(survdiff(Surv(time,status)~group,res.cat)$chisq,1)
  
    pdf(paste0("bestSeparation/",cancer[i],"_",lnc[j],"_survival_sig2.pdf"))
    g<-ggsurvplot(kmfit,risk.table=F,#生存统计统计表
                  # Change legends: title & labels
                  legend.title = "group",
                  legend.labs = c("high","low"),
                  conf.int=TRUE,#添加置信区间带
                  palette = brewer.pal(9,"Set1")[c(3,5)],#颜色设置
                  pval=TRUE,#log-rank检验
                  pval.method=TRUE,#添加检验text
                  pval.size=5,
                  ncensor.plot=FALSE,# 画censor图
                  title=paste(cancer[i],"(",lnc[j],"-",signif(km_p,4),")",sep=''))+
      theme_survminer(font.legend=15,font.x=15,font.y=15,font.main=18)
    print(g)
    dev.off()
    
    km_p_all <- rbind(km_p_all,
                      cbind(cancer[i],lnc[j],km_p))
  }
}
write.table(km_p_all,"bestSeparation/lncRNA_exp_survival_pvalue.txt",sep='\t',quote = F,row.names = F)
km_p_all_sig <- km_p_all[as.numeric(km_p_all[,3]) < 0.05,]
#===========================lncRNA在不同癌症中不同生存组的表达=================================
setwd("J:/xieyunjin/new_result/immune_R1R2R3/example/BLCA-CESC-HNSC-LUSC/BLCA-CESC-LUSC/bestSeparation")
####BLCA#######
cancer <- "BLCA"
lnc <- "ENSG00000262370.1"
group <- read.table(paste0(cancer,"_",lnc,"_bestSeparation.txt"),sep='\t',header = T,as.is = T)

symbol <- read.table("J:/xieyunjin/lncRNA_synergisty_project/pan_cancer_symbol/symbol.txt",sep='\t',header = T,as.is = T)
lnc_symbol <- symbol[symbol$lncRNA %in% lnc,3]
lnc_symbol2 <- gsub("\\..*","",lnc_symbol)

library(ggplot2)
library(RColorBrewer)

pdf(paste0(cancer,"_",lnc,"_exp_boxplot.pdf"))
ggplot(group,aes(x= factor(group,levels = c("low","high")),y=ENSG00000262370.1,fill=group))+
  geom_jitter(aes(color= group),width = 0.2,size = 1)+
  geom_violin(alpha = 0.6)+
  scale_fill_manual(values = c("high"= brewer.pal(9,"Set1")[3],"low"= brewer.pal(9,"Set1")[5]))+
  scale_color_manual(values = c("high"= brewer.pal(9,"Set1")[3],"low"= brewer.pal(9,"Set1")[5]))+
  labs(x = "",y=lnc_symbol2)+
  theme_classic()+
  coord_flip()
dev.off()


####CESC#######
cancer <- "CESC"
lnc <- "ENSG00000234883.2"
group <- read.table(paste0(cancer,"_",lnc,"_bestSeparation.txt"),sep='\t',header = T,as.is = T)

lnc_symbol <- symbol[symbol$lncRNA %in% lnc,3]
lnc_symbol2 <- gsub("\\..*","",lnc_symbol)

pdf(paste0(cancer,"_",lnc,"_exp_boxplot.pdf"))
ggplot(group,aes(x= factor(group,levels = c("low","high")),y=eval(parse(text = lnc)),fill=group))+
  geom_jitter(aes(color= group),width = 0.2,size = 1)+
  geom_violin(alpha = 0.6)+
  scale_fill_manual(values = c("high"= brewer.pal(9,"Set1")[3],"low"= brewer.pal(9,"Set1")[5]))+
  scale_color_manual(values = c("high"= brewer.pal(9,"Set1")[3],"low"= brewer.pal(9,"Set1")[5]))+
  labs(x = "",y=lnc_symbol2)+
  theme_classic()+
  coord_flip()
dev.off()


####LUSC#######
cancer <- "LUSC"
lnc <- "ENSG00000234883.2"
group <- read.table(paste0(cancer,"_",lnc,"_bestSeparation.txt"),sep='\t',header = T,as.is = T)

lnc_symbol <- symbol[symbol$lncRNA %in% lnc,3]
lnc_symbol2 <- gsub("\\..*","",lnc_symbol)

pdf(paste0(cancer,"_",lnc,"_exp_boxplot.pdf"))
ggplot(group,aes(x= factor(group,levels = c("low","high")),y=eval(parse(text = lnc)),fill=group))+
  geom_jitter(aes(color= group),width = 0.2,size = 1)+
  geom_violin(alpha = 0.6)+
  scale_fill_manual(values = c("high"= brewer.pal(9,"Set1")[3],"low"= brewer.pal(9,"Set1")[5]))+
  scale_color_manual(values = c("high"= brewer.pal(9,"Set1")[3],"low"= brewer.pal(9,"Set1")[5]))+
  labs(x = "",y=lnc_symbol2)+
  theme_classic()+
  coord_flip()
dev.off()