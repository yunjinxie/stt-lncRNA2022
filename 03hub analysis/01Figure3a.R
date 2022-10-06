##计算癌症免疫网络中联通度
setwd("J:/xieyunjin/new_result/17_immune_network")
f<-dir("network")
f2<-f[grep("network",f)]
cancer<-gsub("_immune_network.txt","",f2)
p<-file.path("network",f2)
dir.create("degree")
for(i in 1:length(p)){
	net<-read.table(p[i],sep='\t',header=T,as.is=T)
	net2<-as.matrix(net[,c(1,2)])
	degree<-table(c(net2[,1],net2[,2]))
	degree2<-cbind(names(degree),degree)
	write.table(degree2,paste("degree/",cancer[i],"_degree.txt",sep=''),sep='\t',quote=F,row.names=FALSE,col.names=c("lncRNA","degree"))
}
####===========================选择免疫协同网络和17个免疫癌症网络中的hub节点==========================
#####################immune'hub nodes
dir.create("hub/cancer",decursive=TRUE)
immuLnc<-read.table('J:/xieyunjin/new_result/immune/network/network_data/node_attribute_type.txt',sep='\t',header=T)##大网中节点及其连通度
immuLnc2<-as.matrix(immuLnc[,c(1,2)])
hub<-selectHub(immuLnc2,p=0.1)#threshod 
write.table(hub,'hub/immune_big_hub.txt',sep='\t',quote=FALSE,row.names=FALSE)
##17个癌症hub
f<-dir("degree")
cancer<-gsub("_degree.txt","",f)
ad<-file.path("degree",f)

for(i in 1:length(cancer)){
	lnc<-read.table(ad[i],sep='\t',header=T,as.is=T)##lncRNA degree
	lnc$lncRNA<-apply(as.matrix(lnc[,1,drop=F]),1,function(x) substr(x,1,15))
	hub<-selectHub(lnc,p=0.1)##阈值
	write.table(hub,paste("hub/cancer/",cancer[i],"_hub.txt",sep=''),sep='\t',quote=FALSE,row.names=FALSE)
}
###输入节点和相应的连通度
selectHub<-function(nodeDegree,p){
	sortnode<-order(as.numeric(nodeDegree[,2]),decreasing=TRUE)##按照连通度降序排序
	r=round(p*dim(nodeDegree)[1])
	if(r==0){
		return("p is too small")
	}
	if(as.numeric(nodeDegree[sortnode[r],2])!=1){##如果top最后一个节点连通度不为1
	#if the degree of lncRNA in a line equals next one,then also sign the next as a hub 
		while(nodeDegree[sortnode[r],2]==nodeDegree[sortnode[r+1],2]){
		r<-r+1
		}
	hub<-nodeDegree[sortnode[1:r],1]
	degree<-nodeDegree[sortnode[1:r],2]
	de<-cbind(hub,degree)##取出的hub及其连通度
	}else{
		index<-which(as.numeric(nodeDegree[sortnode[1:r],2])!=1)##如果top最后一个节点连通度为1，那么选择连通度不为1的为hub
		if(length(index)==0){
			return("degrees are 1")
		}
		hub<-nodeDegree[sortnode[index],1]
		degree<-nodeDegree[sortnode[index],2]
		de<-cbind(hub,degree)
	}
	colnames(de)<-c("lncRNA","degree")
	return(de)
}
####===================统计免疫大网中hub lncRNA在哪些癌症中也是hub=========================
bighub<-read.table('hub/immune_big_hub.txt',sep='\t',header=T,as.is=T)
bighub2<-as.matrix(bighub[,1,drop=F])
d<-dir("hub/cancer")
cancer<-gsub("_hub.txt","",d)
path<-file.path("hub/cancer",d)

cabighub<-matrix(0,nrow=length(bighub2),ncol=length(cancer))
colnames(cabighub)<-cancer
rownames(cabighub)<-bighub2

for(i in 1:length(path)){
	data<-read.table(path[i],sep='\t',header=T,as.is=T)
	cancerhub<-rownames(data)<-as.matrix(data[,1])
	inter<-intersect(cancerhub,bighub2)
	cabighub[inter,i]<-data[inter,2]
}
write.table(cabighub,"bubble/data/immune_hub_cancer_matrix.txt",sep='\t',quote=F)
##
res<-c()
for(i in 1:dim(cabighub)[1]){
	index<-which(cabighub[i,]!=0)
	cancer<-colnames(cabighub)[index]
	string<-paste(cancer,collapse=',')
	res<-rbind(res,data.frame(bighub2[i],string,length(cancer)))
}
write.table(res,"bubble/data/immune_hub_cancer.txt",sep='\t',quote=F,row.names=F,col.names=c("lncRNA","cancer","length"))
##
common<-res[res[,3]>1,]
specific<-res[res[,3]==1,]
other<-res[res[,3]==0,]
classRes<-Reduce(rbind,list(common,specific,other))
classRes$type<-rep(c("common","specific","other"),
               c(dim(common)[1],dim(specific)[1],dim(other)[1]))
colnames(classRes)<-c("lncRNA","cancer","length","type")
write.table(classRes,"bubble/data/hub_classify.txt",sep='\t',quote=F,row.names=F)
##================免疫大网中hub在各个癌症免疫网络中连通度/节点个数======================
d<-"J:/xieyunjin/new_result/17_immune_network/degree"
f<-dir(d)
path<-file.path(d,f)
cancer<-gsub("_degree.txt","",f)
bighub<-read.table('hub/immune_big_hub.txt',sep='\t',header=T,as.is=T)
colnames(bighub)<-c("lncRNA","degree")
allLnc<-as.matrix(bighub[,1])
per<-matrix(0,nrow=length(allLnc),ncol=length(cancer))
rownames(per)<-allLnc
colnames(per)<-cancer
count<-c()
for(i in 1:length(path)){
	node<-read.table(path[i],sep='\t',header=T,as.is=T)
	lnc<-substr(as.matrix(node[,1]),1,15)
	rownames(node)<-lnc
	inter<-intersect(allLnc,lnc)
	per[inter,i]<-node[inter,2]/length(lnc)
	count<-c(count,dim(node)[1])
}

names(count)<-cancer
write.table(count,"immune_network_node_count.txt",sep='\t',quote=F,col.names="node_num")

##合并type
classRes<-read.table("bubble/data/hub_classify.txt",sep='\t',header=T,as.is=T)
per2<-as.data.frame(per)
per2$lncRNA<-rownames(per2)
result<-Reduce(function(x,y) merge(x,y,by='lncRNA'),list(per2,bighub,classRes))
write.table(result,"bubble/data/hub_classify_cancer_degree_matrix.txt",sep='\t',row.names=F,quote=F)
##=======================免疫hub lncRNA在哪些癌症中协同==========================
setwd("J:/xieyunjin/new_result/17_immune_network/bubble")
hub<-read.table("data/hub_classify.txt",sep='\t',header=T,as.is=T)
hub2<-as.matrix(hub[,1])
d<-"J:/xieyunjin/new_result/17_immune_network/degree"
f<-dir(d)
cancer<-gsub("_degree.txt","",f)
path<-file.path(d,f)
locate<-c()
for(i in 1:length(cancer)){
	canLnc<-read.table(path[i],sep='\t',header=T,as.is=T)
	lnc<-substr(as.matrix(canLnc[,1]),1,15)##一个癌症中协同lncRNA
	hubInter<-intersect(hub2,lnc)
	if(length(hubInter)==0){
		next
	}
	old<-intersect(hubInter,locate[,1])
	if(length(old)==0){
		locate<-rbind(locate,cbind(hubInter,cancer[i]))
	}else{
		index<-which(locate[,1] %in% old)
		locate[index,2]<-apply(locate[index,2,drop=F],1,function(x) paste(c(x,cancer[i]),collapse=','))
		new<-setdiff(hubInter,old)
		if(length(new)==0){
			next
		}
		locate<-rbind(locate,cbind(new,cancer[i]))
	}
}
locate<-as.data.frame(locate)
colnames(locate)<-c("lncRNA","syn_in_cancer")
locate2<-merge(locate,hub,by='lncRNA')
locate2sort<-locate2[order(locate2$length,decreasing=T),]
write.table(locate2sort,"data/hub_classify_syn_in_cancer.txt",sep='\t',quote=F,row.names=F)

#==================免疫hub lncRNA与功能的关系（riverplot）==================
setwd("J:/xieyunjin/new_result/17_immune_network/bubble")

immune_go <- read.table("J:/xieyunjin/lncRNA_synergisty_project/immune_go/immune_go_term2.txt",sep='\t',header=T,as.is = T)

pan_go <- read.table("J:/xieyunjin/new_result/go_lnc_pairs_count/pan_cancer_lnc_go.txt",sep='\t',header = F,as.is = T)

pan_imm_go <- merge(immune_go,pan_go,by.x= "go_names",by.y="V2")

pairs <- do.call(rbind,lapply(as.matrix(pan_imm_go[,"V3"]),function(x) unlist(strsplit(x,"_"))))
#node <- union(pairs[,1],pairs[,2])
#link <- unique(pairs[,])

pan_imm_go_node <- data.frame(unique(rbind(cbind(pan_imm_go[,"type"],pairs[,1]),
                                           cbind(pan_imm_go[,"type"],pairs[,2]))))
head(pan_imm_go_node)
colnames(pan_imm_go_node) <- c("go","lncRNA")

symbol <- read.table("J:/xieyunjin/new_result/pan_cancer_network/data/symbol.txt",sep='\t',header = T,as.is = T)

pan_imm_go_node2 <- merge(pan_imm_go_node,symbol)
########画riverplot图
library(ggalluvial)
library(plyr)

immune_hub <- read.table("data/hub_classify_cancer_degree_matrix_sort.txt",sep='\t',header=T,as.is = T)
immune_hub <- as.matrix(immune_hub)

pan_imm_go_node_hub <- pan_imm_go_node2[pan_imm_go_node2[,3] %in% immune_hub[,"symbol"],]
pan_imm_go_node_hub$go <- gsub("_"," ",pan_imm_go_node_hub$go)

dataP <- data.frame(frenq = 1,
                    Cohort = rep(c(1:nrow(pan_imm_go_node_hub)),times = 2),
                    x = rep(c("go","lncRNA"),each = nrow(pan_imm_go_node_hub)),
                    stratum = c(pan_imm_go_node_hub[,2],pan_imm_go_node_hub[,3]))


immune_col <- read.table("J:/xieyunjin/new_result/immune_type/color.txt",sep='\t',header = F,as.is = T,comment.char='$')
immune_col <- as.matrix(immune_col)
immune_col2 <- immune_col[,1]
names(immune_col2) <- gsub("_"," ",immune_col[,2])

immune_hub2 <- immune_hub[,c("type","symbol")]

color <- rep("#FA8072",times=nrow(immune_hub2))
color[immune_hub2[,"type"] == "specific"] <- "#607B8B"
color[immune_hub2[,"type"] == "other"] <- "#8B7E66"

names(color) <- immune_hub2[,"symbol"]
color2 <- c(immune_col2,color)

pdf("riverplot_2.pdf")
ggplot(dataP,
       aes(x=factor(x,level=c("go","lncRNA")),y=frenq,stratum=factor(stratum,levels = c(immune_hub2[nrow(immune_hub2):1,"symbol"],unique(pan_imm_go_node_hub$go))),alluvium=Cohort,fill=stratum,label=stratum)) +
  geom_flow(width=1/3) +
  geom_stratum(width=1/3,linetype=1,size=0.5,alpha=1,color="white") +
  geom_text(stat="stratum",size=2) +
  scale_x_discrete(limits=c()) +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank()) +
  scale_fill_manual(values=color2)
dev.off()
