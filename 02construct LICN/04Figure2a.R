#================提取泛癌免疫网络===============
goData<-read.table("J:/xieyunjin/lncRNA_synergisty_project/immune_go/immune_go_term2.txt",sep='\t',header=F,as.is=T)
go<-as.matrix(goData[,1])
add<-"J:/xieyunjin/new_result/go_lnc_pairs_count/"
allAD<-paste(add,"pan_cancer_lnc_go.txt",sep='/')
all<-read.table(allAD,sep='\t',header=F,as.is=T)

s<-all[which(as.matrix(all[,2]) %in% as.matrix(go)),]
pairs<-unlist(strsplit(as.character(s[,3]),"_"))
pairs2<-matrix(substr(pairs,1,15),ncol=2,byrow=T)
result<-cbind(s[,2],pairs2)
network<-unique(pairs2)
colnames(network)<-c("lnc1","lnc2")
write.table(result,"network_data/go_pairs.txt",sep='\t',quote=F,row.names=F,col.names=F)
write.table(network,"network_data/immnue_pairs_network.txt",sep='\t',quote=F,row.names=F)

#================节点属性=====================
###节点大小
degree<-sort(table(c(network[,])))
degree2<-as.data.frame(degree)
colnames(degree2)<-c("lncRNA","degree")
###ID2symbol
symbol <- read.table("J:/xieyunjin/lncRNA_synergisty_project/pan_cancer_symbol/symbol.txt",sep='\t',header = T,as.is = T)
symbol <- as.matrix(symbol)

node_symbol <- symbol[match(degree2[,1],symbol[,2]),3]
node_attr <- cbind(degree2,gsub("\\..*","",node_symbol))
colnames(node_attr)[3] <- "symbol" 
write.table(node_attr,"network_data/node_attribute.txt",sep='\t',quote=F,row.names=F)
#================节点在哪些癌症免疫子网出现==================
dd <- "J:/xieyunjin/new_result/17_immune_network/node"
file <- list.files(dd)
path <- file.path(dd,file)
imm_can <- gsub("_node.txt","",file)

can_matrix <- matrix(0,nrow=nrow(node_attr),ncol=length(imm_can))
colnames(can_matrix) <- imm_can
rownames(can_matrix) <- as.matrix(node_attr[,1])

for(i in 1:length(imm_can)){
  data <- read.table(path[i],sep='\t',header = T,as.is = T)
  data <- as.matrix(data)
  lnc <- substr(data[,1],1,15)
  
  inter <- intersect(lnc,node_attr[,1])
  if(length(inter)!=0){
    can_matrix[inter,i] <- 1
  }
  
}

node_attr2 <- cbind(node_attr,can_matrix)

write.table(node_attr2,"network_data/node_attribute-cancer.txt",sep='\t',quote=F,row.names=F)

#####===============边属性，调控免疫功能============================
pan_go<-read.table("F:/xieyunjin/new_result/go_lnc_pairs_count/pan_cancer_lnc_go.txt",sep='\t',header=F,as.is=T)
colnames(pan_go)<-c("GO","go_names","pairs")
immu_go<-read.table("F:/xieyunjin/lncRNA_synergisty_project/immune_go/immune_go_term2.txt",sep='\t',header=T,as.is=T)

immu_pan_go<-merge(pan_go,immu_go,by="go_names")
data<-immu_pan_go[,c(1,3,5)]
write.table(data,"pan_immune_go.txt",sep='\t',quote=F,row.names=F)

immu_net<-read.table("F:/xieyunjin/new_result/immune_type/pan_immune_go.txt",sep='\t',header=T,as.is=T)
result<-ddply(immu_net,.(pairs),summarize,paste(unique(type),collapse=','))
lnc_pairs<-lapply(result$pairs,function(x) unlist(strsplit(x,"_")))
lnc_pairs2<-matrix(unlist(lnc_pairs),nrow=length(lnc_pairs),byrow=T)
result2<-cbind(substr(lnc_pairs2,1,15),result[,2])
write.table(result2,"immune_network_type.txt",sep='\t',quote=F,row.names=F)
