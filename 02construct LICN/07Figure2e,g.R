############构建17个癌症的免疫调控子网，边属性是调控的免疫功能类别
add2<-"F:/xieyunjin/lncRNA_synergisty_project/"
setwd(add2)
immu_go2 <- read.table("immune_go/immune_go_term2.txt",sep='\t',header=T,as.is=T)
####################
add<-"F:/xieyunjin/project_all/project_all/lncRNA_FDR_data/cancer_fit_model_p_adjust"
p<-list.files(add)
cancer<-gsub("_fit_model_p_adjust.txt","",p)
netAD<-file.path(add,p)

setwd("F:/xieyunjin/new_result/17_immune_network/")
dir.create("network")
setwd("network")
for(i in 1:length(cancer)){
	data<-read.table(netAD[i],sep=',',header=T,as.is=T)
	immu_data<-merge(data,immu_go2,by='GO')
	if(dim(immu_data)[1]==0){
		next
	}
	lnc<-unique(immu_data[,c(2,3)])
	lnc_pairs<-paste(lnc[,1],lnc[,2],sep='-')
	immu_data$pairs<-paste(immu_data[,2],immu_data[,3],sep='-')
	lnc_go_all<-c()
	for(j in 1:length(lnc_pairs)){
		index<-which(immu_data$pairs %in% lnc_pairs[j])
		go<-unique(immu_data[index,]$type)
		lnc_go<-cbind(lnc[j,],paste(go,collapse=','))
		lnc_go_all<-rbind(lnc_go_all,lnc_go)
	}
	colnames(lnc_go_all)<-c("lnc1","lnc2","immune_type")
	write.table(lnc_go_all,paste(cancer[i],"_immune_network.txt",sep=''),sep='\t',quote=F,row.names=F)
	write.table(immu_data,paste(cancer[i],"_immune_data.txt",sep=''),sep='\t',quote=F,row.names=F)
}