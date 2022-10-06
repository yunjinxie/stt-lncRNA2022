##### 对之前hyper_test及fdr<0.01的lncRNA对进行 ####
##### 处理：lncRNA对的靶基因集合的交集进行GO功能注释 ##
##### hyper_fdr_*.txt：lncRNA对 ######################
##### linc_gene*.txt：lncRNA及调控的靶基因 ############
##### g2g: GO.BP功能基因集合 ######################
#统计GO里的基因数为 5729个
#N=5729+相应癌症里的miRNA的靶基因数

setwd("J:/xieyunjin/test-GBM-model")
library(data.table)
preg2g<-readLines("J:/xieyunjin/lncRNA_synergisty_project/g2g")
g2g<-strsplit(preg2g,"\t")#编号	mRNA

prelg<-readLines("GBM_convert_lnc_mRNA.txt")
lg<-strsplit(prelg,"\t")
######## 将lg的list格式转化为matrix，所以要填充NA 
plg_matrix<-t(
						sapply(
						lg,
						function(row,max_length)c(row,rep(NA,max_length-length(row))),
						max(sapply(lg,length))
						)
)
lg_matrix<-plg_matrix[,-1]
rownames(lg_matrix)<-plg_matrix[,1]##lncRNA
#####
#fdr<0.01

phf_05<-fread("lnc-lnc-hyper-test01.txt",header=F,data.table=F)
hf_05<-phf_05[,c(1,3)]
hf_05<-as.matrix(hf_05)
#表达谱里的基因数
temp2<-c()
for(we in 1:ncol(lg_matrix)){
temp1<-lg_matrix[,we]
temp2<-c(temp2,temp1)
}
temp3<-unique(temp2)
temp4<-length(temp3)

N<-5729+temp4 #GO和表达谱里的基因的并集
NROW<-nrow(hf_05)
LEN<-length(g2g)

for(i in 1:NROW){
		########### lncRNA基因集
	hi1<-hf_05[i,1]
	hi2<-hf_05[i,2]
	nhf1<-as.character(hi1)
	nhf3<-as.character(hi2)
	prelnc1<-lg_matrix[nhf1,]
	lnc1<-prelnc1[!is.na(prelnc1)]
	prelnc3<-lg_matrix[nhf3,]
	lnc3<-prelnc3[!is.na(prelnc3)]
	########## lncRNA间的基因交集
	prem1<-match(lnc1,lnc3)##返回lnc1在lnc3的参数
	m1<-prem1[!is.na(prem1)]
	geneset<-lnc3[m1]
	GSLEN<-length(geneset)
	
	#############################
	if(GSLEN!=0){
		result<-matrix(nrow=LEN,ncol=7)
		for(j in 1:LEN){
			######### 超几何参数
			m<-GSLEN
			n<-N-m
			k<-length(g2g[[j]])-1

			prem2<-match(geneset,g2g[[j]][-1])
			m2<-prem2[!is.na(prem2)]
			x<-length(m2)##lncRNA间的基因交集和功能集合交集个数
			######### 求p值
			p<-1-phyper(x,m,n,k)
			######### 结果(lnc1,lnc3,lnc1和lnc3交集数目,功能名，功能基因集数目，功能富集数目，p)
			result[j,1]<-hi1
			result[j,2]<-hi2
			result[j,3]<-m
			################### 此处会把001，002，...改为1,2,...
			result[j,4]<-g2g[[j]][1]
			result[j,5]<-k
			result[j,6]<-x
			result[j,7]<-p
		}
		#################### fdr ###########################
		rii<-as.data.frame(result,stringsAsFactors=FALSE)
		#要求功能富集数目>=3
		ri<-rii[rii[,6]>=3,]
		#Rprof("s")
		pfdr<-p.adjust(ri[,7],"fdr")
		log<-pfdr<0.05
		r_log<-ri[log,] #三者有交集并且P小于0.05的lnc-lnc-GO对
		if(nrow(r_log)!=0){	
			pr05<-cbind(r_log,pfdr[log])
			write.table(pr05,"GBM_function_fdr_01_05_1.txt",append=T,quote=F,sep="\t",col.names=F,row.names=F)
			###### fdr<0.01
			pr01<-pr05[pr05[,8]<0.01,]
			###########################################################################################################改
			if(nrow(pr01)!=0){ 
				write.table(pr01,"GBM_function_fdr_01_01_1.txt",append=T,quote=F,sep="\t",col.names=F,row.names=F)
			}
		}
		
	}
}


























