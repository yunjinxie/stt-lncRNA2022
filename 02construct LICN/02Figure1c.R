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

#===============top功能中免疫功能比例======================
setwd("J:/xieyunjin/new_result/go_lnc_pairs_count")

go_count <- read.table("pan_cancer_lnc_go_count.txt",sep='\t',header=F,as.is = T)
#泛癌网络中调控GO功能的lncRNA协同对个数
go_count_sort <- go_count[order(go_count[,3],decreasing = T),]
go_count_sort <- as.matrix(go_count_sort)

immune_go <- read.table("J:/xieyunjin/lncRNA_synergisty_project/immune_go/immune_go_term2.txt",sep='\t',header=T,as.is = T)
immune_go <- as.matrix(immune_go)

####top间隔5####
#top <- c(5,10,15,18,20,25,30,35,40,45,50)/100
#top <- c(5,10,15,18,20,25,30,35,40,45,50)

top <- c(10,15,20,25,30,35,40)
count <- c()
p_value_all <- c()
for(i in top){
  #go_top <- go_count_sort[1:round(nrow(go_count_sort)*i),,drop=F]
  
  go_top <- go_count_sort[1:i,,drop=F]
  n_all <- nrow(go_top)
  
  top_immune <- intersect(immune_go[,4],go_top[,1])
  n_immune <- length(top_immune)
  
  go_non_top <- go_count_sort[(i+1):nrow(go_count_sort),,drop=F]
  n_non_top <- nrow(go_non_top)
  
  non_top_immune <- intersect(immune_go[,4],go_non_top[,1])
  n_non_top_immune <- length(non_top_immune)
  
  
  df <- data.frame(top=c(n_immune,n_all-n_immune),non_top=c(n_non_top_immune,(n_non_top-n_non_top_immune)))
  p_value <- fisher.test(df,alternative = "two.sided")$p.value
  p_value_all <- c(p_value_all,p_value)
  
  count <- rbind(count,cbind(i,n_all,n_immune))
}

write.table(cbind(top,p_value_all),"top_go_immune/top5_fisher.txt",sep='\t',quote = F,row.names = F,col.names = c("top","fisher_pvalue"))

per <- count[,3]/count[,2]
res <- cbind(count,per)

write.table(res,"top_go_immune/top.txt",sep='\t',row.names = F,quote = F)
#n_no_immune
count2 <- cbind(count,count[,2]-count[,3])
plot_data <- cbind(count2[,1,drop=F],count2[,3]/count2[,2],count2[,4]/count2[,2])
#barplot
library(RColorBrewer)
display.brewer.all()

pdf("top_go_immune/top_go_immune_5(number).pdf")

b<- barplot(t(plot_data[,c(2,3)]),col = c(brewer.pal(8,"Set3")[6],brewer.pal(8,"Set2")[7]),xlab = "Top(number)",ylab = "Percentage",
        cex.axis = 1)
axis(1,b,labels = plot_data[,1])

dev.off()

#====================在每个癌症中top前几功能中免疫功能的占比===================
setwd("J:/xieyunjin/new_result/go_lnc_pairs_count")

immune_go <- read.table("J:/xieyunjin/lncRNA_synergisty_project/immune_go/immune_go_term2.txt",sep='\t',header=T,as.is = T)
immune_go <- as.matrix(immune_go[,1])

#癌症go
d <- "J:/xieyunjin/new_result/go_lnc_pairs_count/cancer_go_names"
file <- list.files(d)
path <- file.path(d,file)

cancer <- gsub("_go_names.txt","",file)

top <- c(10,15,20,25,30,35,40)

count_matrix <- matrix(0,nrow=length(top),ncol=length(cancer))
rownames(count_matrix) <- paste("top",top,sep='')
colnames(count_matrix) <- cancer

for(i in 1:length(top)){
  for(j in 1:length(cancer)){
    data <- read.table(path[j],sep='\t',header = T,as.is = T)
    #lncRNA协同对个数
    data <- as.matrix(data)
    
    if(nrow(data) >=top[i]){
      data_sort <- data[order(as.numeric(data[,3]),decreasing = T),]
      data_top <- data_sort[1:top[i],,drop=F]
      
      inter <- intersect(data_top[,2],immune_go)
      count_matrix[i,j] <- length(inter)/top[i]
    }else if(nrow(data)>top[i-1]){
      data_sort <- data[order(as.numeric(data[,3]),decreasing = T),]
    
      inter <- intersect(data_sort[,2],immune_go)
      count_matrix[i,j] <- length(inter)/nrow(data)
    }
  }
}

write.table(count_matrix,"top_go_immune/top_in_cancer.txt",sep='\t',quote = F)


col2 <- colorRampPalette(c("#ffffff","#fffaf0", "#ffffe0", "#fffacd",
                           "#ffe4c4","#ff8c00", "#ff7f50", "#ff0000"))

pdf("top_go_immune/top_immune_in_cancer.pdf",width = 12,height = 7)
corrplot(count_matrix,is.corr=FALSE,method="pie",
         col=col2(100),bg="#fffffc",outline=T,cl.pos='r',
         cl.lim = c(0,1),cl.length =5,cl.cex = 1,cl.ratio = 0.1 ,cl.align = "r",
         tl.cex=1.5,tl.col='black')
dev.off()
