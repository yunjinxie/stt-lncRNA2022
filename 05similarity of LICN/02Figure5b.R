######寻找癌症特异、共享节点和边
setwd("F:/xieyunjin/new_result/immune_R1R2R3")
e0<-list(c("BRCA","LUAD","KIRC","SKCM"),
         c("HNSC","LUSC","BLCA","CESC"),
         c("COAD","UCEC"),
         c("GBM","LGG"),
         c("KIRC","STAD"),c("LUSC","BLCA","CESC"))
         #c("LIHC","PRAD"))
for(i in 1:length(e0)){
  e1<-sort(e0[[i]])
  sePath<-paste(e1,collapse='-')
  dir.create(paste("example",sePath,sep='/'),recursive=T)
  add<-"F:/xieyunjin/new_result/20_immune_network/network"
  file<-dir(add)
  file2<-file[grep(paste(paste(e1,".*","network",sep=''),collapse='|'),file)]
  fileAD<-file.path(add,file2)
  network<-lapply(fileAD,function(x) read.table(x,sep='\t',header=T,as.is=T))
  node<-lapply(network,function(x) unique(c(x$lnc1,x$lnc2)))
  names(node)<-e1
  union_node<-Reduce(union,node)
  result<-matrix(0,nrow=length(union_node),ncol=length(e1))
  rownames(result)<-union_node
  colnames(result)<-e1
  
  edge<-lapply(network,function(x) paste(x$lnc1,x$lnc2,sep='-'))
  names(edge)<-e1
  union_edge<-Reduce(union,edge)
  result2<-matrix(0,nrow=length(union_edge),ncol=length(e1))
  rownames(result2)<-union_edge
  colnames(result2)<-e1
  for(j in 1:length(node)){
    index<-which(union_node %in% node[[j]])
    result[index,j]<-1
    index2<-which(union_edge %in% edge[[j]])
    result2[index2,j]<-1
  }
  setwd(paste("example",sePath,sep='/'))
  write.table(result,"count_node.txt",sep='\t',quote=F)
  write.table(result2,"count_edge.txt",sep='\t',quote=F)
  node_cancer<-unlist(lapply(union_node,function(x) paste(e1[which(result[x,]!=0)],collapse='-')))
  result_node<-cbind(union_node,node_cancer)
  colnames(result_node)<-c("lnc","cancer")
  write.table(result_node,"node.txt",sep='\t',quote=F,row.names=F)
  
  edge_cancer<-unlist(lapply(union_edge,function(x) paste(e1[which(result2[x,]!=0)],collapse='-')))
  union_edge2<-lapply(union_edge,function(x) unlist(strsplit(x,"-")))
  union_edge3<-matrix(unlist(union_edge2),ncol=2,byrow=T)
  result_edge<-cbind(union_edge3,edge_cancer)
  colnames(result_edge)<-c("lnc1","lnc2","cancer")
  write.table(result_edge,"edge.txt",sep='\t',quote=F,row.names=F)
  setwd("F:/xieyunjin/new_result/immune_R1R2R3")
}

##找到共享节点之间的关系
setwd("F:/xieyunjin/new_result/immune_R1R2R3/example/KIRC-STAD")
edge<-read.table("edge.txt",sep='\t',header=T,as.is=T)
node<-read.table("inter/node.txt",sep='\t',header=T,as.is=T,fill=T)
node2<-as.matrix(node[node$type=='lncRNA',1])

index1<-which(as.matrix(edge[,1]) %in% node2)
index2<-which(as.matrix(edge[,2]) %in% node2)
edge_index<-intersect(index1,index2)
edge2<-edge[edge_index,c(1,2)]
write.table(edge2,"inter/edge.txt",append=T,col.names=F,row.names=F,quote=F,sep='\t')

##匹配symbol
setwd("F:/xieyunjin/new_result/immune_R1R2R3/example/BLCA-CESC-LUSC")
data<-read.table("node.txt",sep='\t',header=T,as.is=T)
data<-as.matrix(data)
##all_symbol
sym_data<-read.table("F:/xieyunjin/new_result/pan_cancer_network/data/symbol.txt",sep='\t',header=F,as.is=T)
sym_data<-as.matrix(sym_data)
##
symbol<-apply(data[,1,drop=F],1,function(x) perSym(x,sym_data))
symbol1<-t(symbol)
write.table(symbol1,"node_symbol_gencode.txt",sep='\t',quote=F,row.names=F,col.names=c("lnc","symbol"))
##共享节点
add<-"F:/xieyunjin/new_result/immune_R1R2R3/example/BLCA-CESC-HNSC-LUSC/BLCA-CESC-LUSC/"
node<-read.table(paste(add,"lnc.txt",sep=''),sep='\t',header=F,as.is=T)
node<-as.matrix(node)
symbol<-apply(node,1,function(x) perSym(x,sym_data))
symbol1<-t(symbol)
write.table(symbol1,paste(add,"lnc_symbol.txt",sep=''),sep='\t',quote=F,row.names=F,col.names=c("lnc","symbol"))
##COAD-UCEC
add<-"F:/xieyunjin/new_result/immune_R1R2R3/example/COAD-UCEC"
node<-read.table(paste(add,"node.txt",sep='/'),sep='\t',header=T,as.is=T)
node2<-as.matrix(node[,1,drop=F])
symbol<-apply(node2,1,function(x) perSym(x,sym_data))
symbol1<-t(symbol)
write.table(symbol1,paste(add,"node_symbol_gencode.txt",sep='/'),sep='\t',quote=F,row.names=F,col.names=c("lnc","symbol"))
##
node<-read.table(paste(add,"inter/lnc.txt",sep='/'),sep='\t',header=F,as.is=T)
symbol<-apply(as.matrix(node),1,function(x) perSym(x,sym_data))
symbol1<-t(symbol)
write.table(symbol1,paste(add,"inter/lnc_symbol_gencode.txt",sep='/'),sep='\t',quote=F,row.names=F,col.names=c("lnc","symbol"))

###
perSym<-function(x,sym_data){
  lnc<-x
  index<-grep(lnc,sym_data[,1])
  sym<-sym_data[index,2]
  res<-cbind(x,sym)
  return(res)
}

###寻找共享节点调控的免疫功能并集
library(plyr)
setwd("F:/xieyunjin/new_result/immune_R1R2R3/example")
e1<-dir(getwd())[-7]
file<-c("edge.txt","node.txt")
add<-"F:/xieyunjin/new_result/20_immune_network/network"
file2<-dir(add)
fd<-file2[grep("data",file2)]
fdAD<-file.path(add,fd)
for(i in 1:length(e1)){
  setwd(e1[i])
  cancer<-unlist(strsplit(e1[i],"-"))
  cancerData<-lapply(cancer,function(x) read.table(fdAD[grep(x,fdAD)],sep='\t',header=T,as.is=T))
  nodeGo<-lapply(cancerData,function(x) combGo(x))
  names(nodeGo)<-cancer
  for(j in 1:length(nodeGo)){
    output<-nodeGo[[j]]
    colnames(output)<-c("lnc","go_names")
    write.table(output,paste(cancer[j],"_go_names.txt",sep=''),sep='\t',quote=F,row.names=F)
  }
  node<-read.table("node.txt",sep='\t',header=T,as.is=T)
  result<-ddply(node,.(cancer),function(x) fineNodeFun(x,nodeGo))
  write.table(result,"combine_node_go.txt",sep='\t',quote=F,row.names=F)
  setwd("F:/xieyunjin/new_result/immune_R1R2R3/example")
}
####寻找
combGo<-function(net){#输入一个癌症文件
  colnames(net)<-NULL
  nodeGo<-unique(rbind(cbind(net[,2],net[,5]),cbind(net[,3],net[,5])))
  return(nodeGo)
  
}

fineNodeFun<-function(node.df,node.go){
  type<-node.df[1,2]
  cancer<-unlist(strsplit(as.character(type),"-"))
  can_all<-c()
  for(i in 1:length(cancer)){
    can<-node.go[[cancer[i]]]
    can_all<-unique(rbind(can_all,can))
    colnames(can_all)<-c("lnc","go_names")
  }
  mer1<-merge(node.df,can_all,all.x=TRUE,by.x='lnc')
  return(mer1)
}