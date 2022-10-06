setwd("J:/xieyunjin/new_result/ConsensusClusterPlus")

cluster<-read.table("SKCM_new/cluster2.txt",sep='\t',header=F,as.is=T)
colnames(cluster) <- c("sample","group")
####表达####
lnc<-read.table("SKCM_new/SKCM_filter_lnc.txt",sep='\t',header=F,as.is = F)
lnc<-as.matrix(lnc)

exp<-read.table("J:/xieyunjin/lncRNA_expression/SKCM_regression_1_lncRNA.txt",sep='\t',header=T,as.is = T)
exp2<-exp[,-1]
rownames(exp2)<-substr(as.matrix(exp[,1]),1,15)
#分组
group1<-cluster[cluster$group==1,1]
group2<-cluster[cluster$group==2,1]

lnc_exp1<-exp2[lnc,group1]
lnc_exp2<-exp2[lnc,group2]

mean1<-rowMeans(lnc_exp1)
mean2<-rowMeans(lnc_exp2)

mm <- cbind(mean1,mean2)

pvalue<-t.test(mean1,mean2)$p.value
write.table(pvalue,"SKCM_new/SKCM_two_subtype_expression_t.test.txt",sep='\t',quote=F,col.names = F,row.names = F)
dataP<-data.frame(exp=c(mean1,mean2),
                  group=rep(c("group1","group2"),c(length(mean1),length(mean2))))

cols<-brewer.pal(9,"Set1")[c(3,5)]
names(cols)<-c("group1","group2")

pdf("SKCM_new/SKCM_two_subtype_expression.pdf",width = 7,height = 7)
ggplot(dataP,aes(x=group,y=exp,fill=group))+geom_violin()+geom_boxplot(width=0.2,fill="white")+
  geom_jitter(size=2.5,position=position_jitter(0.2),color="#8a2be2")+
  scale_fill_manual(values=cols)+
  labs(x='',y='Expression(mean) of lncRNAs')+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=12,color="black"))+
  theme(axis.text.x = element_text(size=15,color='black'))+
  theme(axis.title.y = element_text(size=15,color='black'))+
  theme(legend.position = 'none')+
  geom_signif(comparisons = list(c("group1","group2")),test="t.test",tip_length = 0,
              map_signif_level = T,textsize = 6)+
  theme(plot.margin = unit(c(1,2,1,2),"cm"))

dev.off()  
