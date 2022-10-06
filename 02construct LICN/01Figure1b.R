#=============================泛癌网络======================================================
setwd("F:/xieyunjin/new_result/pan_cancer_network")
data<-read.table("bignetwork.txt",sep='\t',header=T)
data<-as.matrix(data)
#==============================节点出现在几种癌症中=========================================
node<-unique(as.vector(data[,c(1,2)]))
d2<-"D:/xieyunjin/project_all/project_all/zhibiao/计算d值/"
can<-read.table(paste(d2,"compute1.txt",sep=''),sep='\t',header=T)
can2<-can[,-1]
rownames(can2)<-can[,1]
title2<-cbind("lncRNA","degree","cancer_length","cancer")
write.table(title2,"pan_cancer_network_node.txt",sep='\t',row.names=F,col.names=F,quote=F,append=T)

degree<-table(as.vector(data[,c(1,2)]))
for(i in 1:length(node)){
   cancer<-colnames(can2)[which(can2[node[i],]!=0)]
   can_len<-length(cancer)
   result<-cbind(node[i],degree[node[i]],can_len,paste(cancer,collapse=','))
   write.table(result,"pan_cancer_network_node.txt",sep='\t',row.names=F,col.names=F,quote=F,append=T)

}
#===================================================ID2symbol=============================
symbol<-read.table("F:/xieyunjin/协同调控lncRNA课题/lncipedia_5_2_hg38_symbol_id.txt",sep='\t',header=T)

s1<-c()
for(i in 1:length(node)){
   
      index<-which(symbol[,2] %in% substr(node[i],1,15))
      if(length(index)!=0){
          if(length(index)==1){
             s1<-rbind(s1,cbind(node[i],as.matrix(symbol[index,1])))
          }else{
             s1<-rbind(s1,cbind(node[i],as.matrix(symbol[index[1],1])))
          }
          
      }else {
         s1<-rbind(s1,cbind(node[i],node[i]))
      }
}
colnames(s1)<-c("lncRNA","symbol")

#=============================================piechart=======================================
cancer_color<-c(BLCA='#8057B5',BRCA='#57B5E5',CESC='#EFB1BE',COAD='#E41E25',GBM='#FDBF6D',HNSC='#B3D789',KICH='#E967F4',KIRC='#F59999',KIRP='#F57F21',LGG='#CAB3D5',LIHC='#A4CEE4',LUAD='#B15927',LUSC='#D3EDEF',OV='#F9F385',PRAD='#F79CD8',READ='#93BA90',SKCM='#E2DB59',STAD='#6EB6BC',THCA='#B25BA8',UCEC='#30A147')

data<-read.table("pan_cancer_network_node.txt",sep='\t',header=T)
data<-as.matrix(data)

attribute_all<-c()
value_matrix<-matrix(0,nrow=dim(data)[1],ncol=length(cancer_color))
colnames(value_matrix)<-names(cancer_color)
for(i in 1:dim(data)[1]){
   cancer<-as.character(unlist(strsplit(data[i,4],",")))
   
   value_matrix[i,cancer]<-100/length(cancer)
   attribute<-paste("piechart:attributelist=",'"',data[i,4],'"',' colorlist="',
                    paste(cancer_color[cancer],collapse=',',sep=''),'"'," showlabels=false",sep='')
   
   attribute_all<-rbind(attribute_all,attribute)
}


data2<-as.data.frame(data)
data3<-cbind(data2,value_matrix)
data3$Gradient<-attribute_all
data3$symbol<-s1[,2]
num<-as.matrix(data3$cancer_len)
data3$symbol2<-data3$symbol
data3$symbol2[as.numeric(num)<4]<-rep("",length(which(as.numeric(num)<4)))
write.table(data3,"pan_cancer_network_node_attribute2.txt",quote=F,sep='\t',row.names=F)
#==================================multi color render=============================
data<-read.table("pan_cancer_network_node_attribute2copy.txt",sep='\t',header=T,strip.white = TRUE)
data<-as.matrix(data)
all_color<-read.table("F:/xieyunjin/lncRNA_synergisty_project/20_color.txt",sep='\t',header=F,comment.char = "/")
all_color<-as.matrix(all_color)
colorList<-c()
for(i in 1:dim(data)[1]){
   cancer<-as.matrix(unlist(strsplit(data[i,4],",")))
   colorList<-c(colorList,
               paste(all_color[which(all_color[,1] %in% cancer),3],collapse=";"))
   
}
result<-cbind(data[,1],colorList)
colnames(result)<-c("Id","colourList")
write.table(result,"pan_cancer_network_node_colorList.txt",sep='\t',quote=F,row.names=F)