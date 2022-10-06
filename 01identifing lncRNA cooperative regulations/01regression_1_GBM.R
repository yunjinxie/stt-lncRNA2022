#regression_1:gene-methylation-cnv-lncRNA
#GBM
GBM_gene<-read.table("/pub4/xujuan/wgj/project2015/results/regression_1_mRNA_data/GBM_regression_1_mRNA.txt",header=T,sep="\t")
GBM_methylation<-read.table("/pub4/xujuan/wgj/project2015/results/regression_1_methylation_data/GBM_regression_1_methylation.txt",header=T,sep="\t")
GBM_cnv<-read.table("/pub4/xujuan/wgj/project2015/results/regression_1_CN_data/GBM_regression_1_CN.txt",header=T,sep="\t")

GBM_lncRNA<-read.table("/pub4/xujuan/wgj/project2015/results/regression_1_lncRNA_data/filter_lncRNA_log2/GBM_regression_1_lncRNA_end.txt",header=T,sep="\t")


methylation_t1<-GBM_methylation[,-1]
methylation_t2<-which(apply(methylation_t1,1,function(x)all(is.na(x)))) # the lines
methylation_t3<-GBM_methylation[-methylation_t2,]
dim(methylation_t3) #20083
dim(GBM_methylation) #21052



#log2  gene 
lo_1<-GBM_gene[,-1]
lo_2<-log(lo_1,2)
lo_3<-GBM_gene[,1]
lo_3<-as.matrix(lo_3) # gene_id
colnames(lo_3)<-"gene_symbol"
lo_4<-cbind(lo_3,lo_2)
GBM_gene_1<-lo_4




gene_id_gene<-GBM_gene_1[,1]
gene_id_meth<-methylation_t3[,1]
gene_id_cnv<-GBM_cnv[,1]

gene_id_1<-intersect(gene_id_gene,gene_id_meth)
#13514
gene_id_2<-intersect(gene_id_1,gene_id_cnv)
#13253

gene_id<-as.matrix(gene_id_2)

colnames(gene_id)<-"gene_id"
# the same gene:gene vs cnv

cnv_1<-merge(gene_id,GBM_cnv,by="gene_id")

# the same gene:gene vs methylation
colnames(gene_id)<-"gene_symbol"
methylation_1<-merge(gene_id,GBM_methylation,by="gene_symbol")
gene_1<-merge(gene_id,GBM_gene_1,by="gene_symbol")


dim(gene_1) 
dim(cnv_1)
dim(methylation_1)
#12770 250 
# the gene number is 12770

# the rownames(gene) is consistant: gene,methylation,cnv

gene_1<-as.matrix(gene_1)
methylation_1<-as.matrix(methylation_1)
cnv_1<-as.matrix(cnv_1)
GBM_lncRNA<-as.matrix(GBM_lncRNA)


all_result_methylation<-c()
all_result_cnv<-c()
all_result_lncRNA<-c()
all_result_intercept<-c()


NROWG<-nrow(gene_1)  # the number of genes
NROWL<-nrow(GBM_lncRNA)  # the number of lncRNAs

# split in five part to run it
##############################
#1:2000 genes
for(j in 1:2000){
    result<-c()
	l_gene<-gene_1[j,-1]
	l_methylation<-methylation_1[j,-1]
	l_cnv<-cnv_1[j,-1]
		
	l_gene<-as.numeric(l_gene)
	gene_inf_1<-which(l_gene=="-Inf")
    l_gene[gene_inf_1]<-10^(-32)
    l_methylation<-as.numeric(l_methylation)
    l_cnv<-as.numeric(l_cnv)	
	for(i in 1:NROWL){

    l_lnc<-GBM_lncRNA[i,-1]
	l_lnc<-as.numeric(l_lnc)
	lnc_inf_1<-which(l_lnc=="-Inf")
    l_lnc[lnc_inf_1]<-10^(-32)
	l<-lm(l_gene~l_methylation+l_cnv+l_lnc,data.frame(l_gene,l_methylation,l_cnv,l_lnc))
	s<-summary(l)
	
	result_s_l<-matrix(nrow=1,ncol=10)
	
	result_s_l[1,1]<-gene_id[j]
	result_s_l[1,2]<-GBM_lncRNA[i,1]
	result_s_l[1,3]<-s$coefficients[2] #methylation coefficients
	result_s_l[1,4]<-s$coefficients[14] #methylation P value
	result_s_l[1,5]<-s$coefficients[3]  #CNV coefficients
	result_s_l[1,6]<-s$coefficients[15]  #CNV P value
	result_s_l[1,7]<-s$coefficients[4] #lncRNA coefficients
	result_s_l[1,8]<-s$coefficients[16] #lncRNA P value
	result_s_l[1,9]<-s$coefficients[1]  # Intercept 
	result_s_l[1,10]<-s$coefficients[13]  #Intercept P value
	result<-rbind(result,result_s_l)
	}

	
	
	###########（gene_id lnc_id methylation_coefficients methylation P value CNV coefficients CNV P value lncRNA coefficients lncRNA P value  Intercept 	  Intercept P value）
	
	# fdr p value 
	gene_lnc_meth_1<-result[,c(1:4)]
	gene_lnc_cnv_1<-result[,c(1,2,5,6)]
	gene_lnc_lnc_1<-result[,c(1,2,7,8)]
	gene_lnc_intercept_1<-result[,c(1,2,9,10)]
	
	p_meth<-gene_lnc_meth_1[,4]
	p_cnv<-gene_lnc_cnv_1[,4]
	p_lnc<-gene_lnc_lnc_1[,4]
	p_intercept<-gene_lnc_intercept_1[,4]
	
	fdr_meth<-p.adjust(p_meth,method="fdr")
	fdr_cnv<-p.adjust(p_cnv,method="fdr")
	fdr_lnc<-p.adjust(p_lnc,method="fdr")
	fdr_intercept<-p.adjust(p_intercept,method="fdr")
	
    fdr_meth<-as.matrix(fdr_meth)
    fdr_cnv<-as.matrix(fdr_cnv)
	fdr_lnc<-as.matrix(fdr_lnc)
	fdr_intercept<-as.matrix(fdr_intercept)
	gene_lnc_meth_2<-cbind(gene_lnc_meth_1,fdr_meth)
	gene_lnc_cnv_2<-cbind(gene_lnc_cnv_1,fdr_cnv)
	gene_lnc_lnc_2<-cbind(gene_lnc_lnc_1,fdr_lnc)
	gene_lnc_intercept_2<-cbind(gene_lnc_intercept_1,fdr_intercept)
	colnames(gene_lnc_meth_2)<-c("gene_id","lncRNA_id","methylation_coefficients","methylation_Pvalue","methylation_FDR")
	colnames(gene_lnc_cnv_2)<-c("gene_id","lncRNA_id","cnv_coefficients","cnv_Pvalue","cnv_FDR")
	colnames(gene_lnc_lnc_2)<-c("gene_id","lncRNA_id","lncRNA_coefficients","lncRNA_Pvalue","lncRNA_FDR")
	colnames(gene_lnc_intercept_2)<-c("gene_id","lncRNA_id","intercept_coefficients","intercept_Pvalue","intercept_FDR")
		
	
	
	 # output the correclation which the fdr<0.1.
	k_meth<-which(fdr_meth<0.1)
	k_cnv<-which(fdr_cnv<0.1)
	k_lnc<-which(fdr_lnc<0.1)
	k_intercept<-which(fdr_intercept<0.1)
	
	result_methylation<-gene_lnc_meth_2[k_meth,]

	result_cnv<-gene_lnc_cnv_2[k_cnv,]

	result_lncRNA<-gene_lnc_lnc_2[k_lnc,]
	
	result_intercept<-gene_lnc_intercept_2[k_intercept,]
	
	
	all_result_methylation<-rbind(all_result_methylation,result_methylation)
	all_result_cnv<-rbind(all_result_cnv,result_cnv)
	all_result_lncRNA<-rbind(all_result_lncRNA,result_lncRNA)
	all_result_intercept<-rbind(all_result_intercept,result_intercept)
}


write.table(all_result_methylation,"/pub4/xujuan/wgj/project2015/results/regression_1_result/GBM/GBM_lm_1_methylation_result_part1.txt",row.names=F,sep="\t",quote=F)	

write.table(all_result_cnv,"/pub4/xujuan/wgj/project2015/results/regression_1_result/GBM/GBM_lm_1_cnv_result_part1.txt",row.names=F,sep="\t",quote=F)	

write.table(all_result_lncRNA,"/pub4/xujuan/wgj/project2015/results/regression_1_result/GBM/GBM_lm_1_lncRNA_result_part1.txt",row.names=F,sep="\t",quote=F)	

write.table(all_result_intercept,"/pub4/xujuan/wgj/project2015/results/regression_1_result/GBM/GBM_lm_1_intercept_result_part1.txt",row.names=F,sep="\t",quote=F)	
