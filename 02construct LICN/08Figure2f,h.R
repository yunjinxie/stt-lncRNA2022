###==================调控免疫功能lncRNA协同对个数====================
library(data.table)

setwd("F:/xieyunjin/new_result/go_lnc_pairs_count")
####################17个癌症网络中所有调控的免疫功能
add<-"F:/xieyunjin/project_all/project_all/lncRNA_FDR_data/cancer_fit_model_p_adjust"
p<-dir(add)
pAD<-paste(add,p,sep='/')
cancer<-gsub("_fit_model_p_adjust.txt","",p)
all_go<-c()
dir.create("cancer_immune_count")

immune<-read.table("F:/xieyunjin/lncRNA_synergisty_project/immune_go/immune_go_term2.txt",sep='\t',header=T,as.is=T)
immune2<-immune[order(immune$type),]
id2f<-read.table("F:/xieyunjin/lncRNA_synergisty_project/id2function_GO_Id",sep='\t',header=F,as.is=TRUE)
colnames(id2f)<-c("go_names","GO","go_terms")

for(i in 1:length(cancer)){
	can<-fread(pAD[i],sep=",",header=T)#lnc1,lnc2,GO,p_adjust
   go<-unique(can$GO)
   immu_go<-intersect(as.matrix(go),as.matrix(immune2[,4]))
	index<-which(immune2$GO %in% immu_go)
	immu_names<-immune2[index,1]
	count<-can[,dim(.SD)[1],by=GO]
	colnames(count)<-c("GO","lncRNA_num")##lncRNA协同对个数
    count2<-merge(id2f[,c(1,2)],count,by="GO")#3列，go_names,GO,lncRNA_num
	rownames(count2)<-count2$go_names
	count3<-count2[immu_names,]	
   write.table(count3,paste("cancer_immune_count/",cancer[i],"_lnc_count.txt",sep=''),sep='\t',quote=F,row.names=F)
}

###################
#========================================画柱状图============================================
library(ggplot2)
add2<-"cancer_immune_count"
p<-dir(add2)
pAD<-paste(add2,p,sep='/')#分别统计的17个癌症中调控GO功能的lncRNA协同对个数
###############################################
cancer<-gsub("_lnc_count.txt","",p)
can_col<-read.table("F:/xieyunjin/lncRNA_synergisty_project/20_color.txt",sep='\t',header=F,comment.char="/",as.is=TRUE)
dir.create("immune_result_type_new")
output<-paste("immune_result_type_new/",cancer,"_go_barplot.pdf",sep='')
par(mai=c(1,1,1,1))
for(i in 1:length(cancer)){
	dataP<-read.table(paste("cancer_immune_count/",cancer[i],"_lnc_count.txt",sep=''),sep='\t',header=T,as.is=T)
	dataP$go_names<-tolower(gsub("_"," ",dataP$go_names))
	##ggplot
    pdf(output[i],width=10)
	
    p<-ggplot(dataP,aes(x=factor(go_names,levels=go_names),y=lncRNA_num))+geom_bar(stat="identity",fill=can_col[i,2],col="gray")+
    labs(y="number(pairs)",x="")+
    coord_flip()+
      theme_bw()+
    theme(axis.text.y=element_text(size=if(dim(dataP)[1]<50){11} else if(dim(dataP)[1]<90){10} else{9},debug=T,hjust=.5))
      
    print(p)
   
    dev.off()
}
