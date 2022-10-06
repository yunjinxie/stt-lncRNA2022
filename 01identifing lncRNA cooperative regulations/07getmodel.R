### ���������繦��ģ�� 
### ����ģ��������㣺
### 1.ģ���а���>=3������ 
### 2.���������������ʣ�
### 	1>ģ���������������<=2 
### 	2>����·������CPL С���漴���� 
### ע�⣺���������粢���������������нڵ�
#################################################
#��������pֵ��û�н���FDRУ��
#################################################
### GO���ܼ���#���ܸ���fdr<0.05

setwd("J:/xieyunjin/test-GBM-model")
library(data.table)

gene_id<-fread("J:/xieyunjin/test-GBM-model/data/gene_id_symbol.txt",sep=",",data.table=F)
### GO���ܼ���
pregene2go<-readLines("J:/xieyunjin/lncRNA_synergisty_project/g2g")
g2g<- strsplit(pregene2go, "\t")
g2g_matrix <- t(
						sapply(
						g2g,
						function(row, max_length)c(row, rep(NA, max_length - length(row))),
						max(sapply(g2g, length))
						)
)
rownames(g2g_matrix) <- g2g_matrix[, 1]
### lncRNA���صĻ��򼯺� 
prelg<-readLines("GBM_convert_lnc_mRNA.txt")
lg <- strsplit(prelg, "\t")
### ��lg��list��ʽת��Ϊmatrix������Ҫ���NA 
lg_matrix <- t(
						sapply(
						lg,
						function(row, max_length)c(row, rep(NA, max_length - length(row))),
						max(sapply(lg, length))
						)
)
rownames(lg_matrix) <- lg_matrix[, 1]
lg_matrix=lg_matrix[,-1]
### ����GO���ܸ�������ļ� 
lnci_lncj_func <- fread("GBM_function_fdr_01_01_1.txt", header=F,data.table=F)
lnci_lncj_func<-as.matrix(lnci_lncj_func)
################################################
### ���� ########################################
library(igraph)
### ���뵰�������� 
mcp <- fread("J:/xieyunjin/test-GBM-model/data/maxcomponent7.28.txt", header=F,data.table=F)
PPI <- graph.data.frame(mcp, directed=F)
### �����������,��ת����igraphͼ����
random_net_names <- dir("J:/xieyunjin/test-GBM-model/data/randomnetwork")
random_nets <- lapply(random_net_names, function(x){
												path <- paste("J:/xieyunjin/test-GBM-model/data/randomnetwork", x, sep="/")
												rn <- fread(path, header = F,data.table=F)
												graph.data.frame(rn, directed = F)
										}
					)
### ���綥����
vertice_names <- get.data.frame(PPI, what = "vertices")[,1]
################################################
### �����򣬼���ģ���� ############################
NROW <- nrow(lnci_lncj_func)
### ������ļ�
outfile <- file("GBM_fit_model_p_2.txt", "w")
for(each_row in 1 : NROW){
	row_data <- lnci_lncj_func[each_row, ]
	### row_data�ĵ�6��(ģ���ѡ������) �� ��3�У�lncRNA�Թ�ͬ�л���ĸ������ı��� >= 0.01
	opp1<-as.numeric(row_data[6])
	opp2<-as.numeric(row_data[3])
	if((opp1 / opp2) >= 0.01){
		input <- as.character(row_data[c(1, 2, 4)])
		### lncRNA_i ���򼯺� 
		lnci <- input[1]
		lnci_genes <- lg_matrix[lnci, ]
		i_g <- lnci_genes[! is.na(lnci_genes)]
		### lncRNA_j ���򼯺�
		lncj <- input[2]
		lncj_genes <- lg_matrix[lncj, ]
		j_g <- lncj_genes[! is.na(lncj_genes)]
		### ��ѡģ�� ���򼯺�
		fc <- input[3]
		fc<-as.numeric(fc)
		fc<-as.character(fc)
		lncij_genes <- intersect(i_g, j_g)
		Gs_id <- intersect(g2g_matrix[fc, - 1 ], lncij_genes)
		#��Gs_id��symbolת��entrze id
		temp_G1<-as.matrix(Gs_id)
        colnames(temp_G1)<-"HGNC.symbol"
        temp_G2<-merge(temp_G1,gene_id,by="HGNC.symbol")
        temp_G3<-temp_G2[,2]
		Gs<-unique(temp_G3)
		Gs<-as.character(Gs)
		### ���� Gs�Ƿ�ȫ������������
		test_Gs <- all(! is.na(match(Gs, vertice_names)))##all()����Ƿ�ȫ��Ϊ�棬��ȫΪ�淵��TURE
		if(test_Gs){
			##############################################
			### ���� �����ʻ������� ģ�������������ľ���D1�Ƿ�<3
			### �Լ��� �����ʻ������� �е�cpl D2�Ƿ�С���������
			D1 <- shortest.paths(PPI, Gs, Gs)
			if(all(D1 < 3)){
				len_Gs <- length(Gs) 
				num_pairs <- choose(len_Gs, 2) * 2
				D2 <- sum(D1) / num_pairs 
				### ������򼯺�Gs��������������CPL
				CPLs <- lapply(random_nets, function(r_net){
														sum(shortest.paths(r_net, Gs, Gs)) / num_pairs
														}
								       )
				### ���� pֵ
				p_value <- sum(CPLs < D2) / 1000
				out_value<-c(input, p_value) 
				cat(out_value, sep = "\t", file = outfile)
					cat("\n", file = outfile)
			}
		}
	}
}
close(outfile)

















