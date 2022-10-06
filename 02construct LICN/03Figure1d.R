#================================统计20个癌症中调控不同go功能的lncRNA协同对个数=======================
library(data.table)

setwd("J:/xieyunjin/new_result/go_lnc_pairs_count")
####################20个癌症网络中所有调控的功能
add<-"F:/xieyunjin/project_all/project_all/lncRNA_FDR_data/cancer_fit_model_p_adjust"
p<-dir(add)
pAD<-paste(add,p,sep='/')
cancer<-gsub("_fit_model_p_adjust.txt","",p)
all_go<-c()
dir.create("cancer_count")
for(i in 1:length(cancer)){
	can<-fread(pAD[i],sep=",",header=T)#lnc1,lnc2,GO,p_adjust
   go<-unique(can$GO)
	count<-can[,dim(.SD)[1],by=GO]
   colnames(count)<-c("GO","lncRNA_num")##lncRNA协同对个数
   write.table(count,paste("cancer_count/",cancer[i],"_lnc_count.txt",sep=''),sep='\t',quote=F,row.names=F)
	all_go<-unique(c(all_go,go))
}
write.table(all_go,"all_go.txt",quote=F,sep='\t',col.names=F,row.names=F)
########
add<-"F:/xieyunjin/new_result/hub_go"
id2f<-read.table(paste(add,"id2function_GO_Id",sep='/'),sep='\t',header=F,as.is=TRUE)
#===========================不分癌症统计泛癌网络中调控GO功能的lncRNA协同对个数============================
all_go<-read.table("all_go.txt",sep='\t',header=F,as.is=TRUE)
colnames(all_go)<-"GO"
all_go_names<-merge(id2f[,c(1,2)],all_go,by="GO")#两列，go_names,GO
add3<-"F:/xieyunjin/project_all/project_all/lncRNA_FDR_data/cancer_fit_model_p_adjust"
p<-dir(add3)
pAD<-paste(add3,p,sep='/')
cancer<-gsub("_fit_model_p_adjust.txt","",p)
result_all<-c()
for(i in 1:length(cancer)){
	can_go<-fread(pAD[i],sep=',',header=T)
	setkey(can_go,"GO")
	
	for(j in 1:dim(all_go_names)[1]){
		go<-all_go_names[j,1]
		go_lnc<-as.matrix(can_go[.(go),.(lnc1,lnc2),nomatch=0])
		pairs<-paste(go_lnc[,1],go_lnc[,2],sep='_')
		if(i==1){
			result_all<-c(result_all,list(pairs))
		}else{
			old<-result_all[[j]]
			inter<-intersect(pairs,old)
			if(length(inter)==0){
				new<-c(old,pairs)
				result_all[[j]]<-new
			}else{
				new<-setdiff(pairs,inter)
				all<-c(old,new)
				result_all[[j]]<-all
			}
		}
	}
}
for(k in 1:length(result_all)){
	count<-length(result_all[[k]])
    go<-all_go_names[k,1]
    go_names<-all_go_names[k,2]
	output<-cbind(rep(go,count),rep(go_names,count),as.matrix(result_all[[k]]))
	output2<-cbind(all_go_names[k,],count)
	write.table(output,"pan_cancer_lnc_go.txt",sep='\t',append=T,row.names=F,col.names=F,quote=F)
	write.table(output2,"pan_cancer_lnc_go_count.txt",quote=F,sep='\t',append=T,row.names=F,col.names=F)

}

#==============柱状图===========
###挑选被最多lncRNA调控的18个GO功能展示
setwd("J:/xieyunjin/new_result/go_lnc_pairs_count")

library(RColorBrewer)
display.brewer.all()

dataP<-read.table("pan_cancer_lnc_go_count_plot.txt",sep='\t',header=T,as.is=T,colClass=c("character","numeric","character"))
# ggplot(dataP,aes(x=go_names,y=count,fill=immune))+geom_bar(stat="identity")+
# coord_flip()+
# theme(axis.text.y=element_text(size=10,debug=T,hjust=.5))+
# scale_fill_manual(values=c("yes"=))
####barplot
data<-as.matrix(dataP[,2])
rownames(data)<-dataP[,1]
plot(0,ylim=c(1,dim(dataP)[1]),xlim=c())
len<-dim(dataP)[1]

dataP$color[dataP$immune=="yes"]<-brewer.pal(8,"Set3")[6]
dataP$color[dataP$immune=="no"]<-brewer.pal(8,"Set2")[7]

cols<-as.vector(dataP$color)

# dataP$density[dataP$immune=="yes"]<-40
# dataP$density[dataP$immune=="no"]<--1
# density<-as.vector(dataP$density)

dataP[,1] <- gsub("_"," ",dataP[,1])

pdf("go_bubble/pan_cancer_go_new.pdf")
par(mai=c(1,3.5,1,.1))
b<-barplot(data,horiz=T,col=cols,beside=T,
           #density=density,
           las=2,xaxt='n',names.arg=dataP[,1])
axis(1,las=1)
legend(200,22,xpd=T,legend=c("immune","non_immune"),fill=cols,
       #density=c(40,-1),
       bty='n')
dev.off()
####################par恢复默认
opar<-par(no.readonly = TRUE)
par(opar)

###============================泛癌画气泡图，行是功能，列是癌症=============
library(data.table)
library(plyr)
setwd("J:/xieyunjin/new_result/go_lnc_pairs_count/go_bubble")
###pan_cancer_lnc_go_count_plot.txt将go name对应到id2function_GO_Id上的go编号后的文件
go<-read.table("plot_go.txt",sep="\t",header=T,as.is=T,colClass=c("numeric","character"))
add<-"J:/xieyunjin/project_all/project_all/lncRNA_FDR_data/cancer_fit_model_p_adjust"
files<-list.files(add)
canAD<-paste(add,files,sep='/')
cancer<-gsub("_fit_model_p_adjust.txt","",files)
plot_matrix<-matrix(0,nrow=dim(go)[1],ncol=length(cancer))
rownames(plot_matrix)<-go[,2]
colnames(plot_matrix)<-cancer
title<-cbind("lnc1","lnc2","GO","cancer")
write.table(title,"result.txt",append=T,quote=F,row.names = F,col.names=F,sep='\t')
for(i in 1:length(cancer)){
	can<-fread(canAD[i],sep=',',header=T)
	setkey(can,"GO")
	can2<-can[.(go[,1]),-4,nomatch=0]
	if(nrow(can2)!=0){
	  output<-cbind(can2,cancer[i])
	  write.table(output,"result.txt",append=T,quote=F,row.names = F,col.names=F,sep='\t')
	  count<-can2[,dim(.SD)[1],by="GO"]
	  colnames(count)<-c("GO","pairs_num")
	  count2<-merge(go,count,by="GO")
	  plot_matrix[as.matrix(count2[,2]),i]<-count2[,3]
	}
}
write.table(plot_matrix,"go_cancer_matrix.txt",sep='\t',quote=F)
##
total<-read.table("result.txt",sep='\t',header=T,as.is=T)
total$pairs<-paste(total[,1],total[,2],sep='-')
count<-ddply(total,.(GO),function(x) length(unique(x$pairs)))

#===============bubble===========================
setwd("J:/xieyunjin/new_result/go_lnc_pairs_count/go_bubble")
library(ggplot2)
library(RColorBrewer)
plot_matrix<-read.table("go_cancer_matrix.txt",sep='\t',header=T,as.is=T)
plot_matrix<-as.matrix(plot_matrix)
rs<-colSums(plot_matrix)
#dataP<-plot_matrix[,rs!=0]
dataP<-plot_matrix
rownames(dataP)<-tolower(gsub("_"," ",rownames(dataP)))
display.brewer.all()
rlen<-dim(dataP)[1]
clen<-dim(dataP)[2]
x1<-rep(1:clen,each=rlen)
y1<-rep(1:rlen,times=clen)
r<-as.numeric(dataP)
r[r==0]<-NA
dataP2<-data.frame(x=x1,y=y1,num=r)
cancer<-colnames(dataP)
##
pdf("go_cancer_bubble2_new.pdf",width = 10,height=5)

ggplot(dataP2,aes(x,y))+geom_point(aes(size=num,col=num))+
  theme_bw()+
  scale_x_continuous(breaks=c(1:length(cancer)),labels=cancer)+
   scale_y_continuous(breaks=c(1:rlen),limits=c(1,rlen),labels=rownames(dataP))+labs(x="",y="")+
 theme(axis.text.y=element_text(size=rel(1.5),hjust = 1))+
 theme(axis.ticks.length =unit(0,"cm"))+
 scale_colour_gradient(low=brewer.pal(9,"OrRd")[3],high=brewer.pal(9,"OrRd")[9])+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        panel.grid = element_blank())
dev.off()

