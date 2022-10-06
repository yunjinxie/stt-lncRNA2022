### 蛋白质网络功能模块 
### 功能模块必需满足：
### 1.模块中包含>=3个基因 
### 2.满足两个拓扑性质：
### 	1>模块内两两基因距离<=2 
### 	2>特征路径长度CPL 小于随即网络 
### 注意：蛋白质网络并不包括所检测的所有节点
#################################################
#输出结果的p值，没有进行FDR校正
#################################################
### GO功能集合#功能富集fdr<0.05

setwd("J:/xieyunjin/test-GBM-model")
library(data.table)

gene_id<-fread("J:/xieyunjin/test-GBM-model/data/gene_id_symbol.txt",sep=",",data.table=F)
### GO功能集合
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
### lncRNA调控的基因集合 
prelg<-readLines("GBM_convert_lnc_mRNA.txt")
lg <- strsplit(prelg, "\t")
### 将lg的list格式转化为matrix，所以要填充NA 
lg_matrix <- t(
						sapply(
						lg,
						function(row, max_length)c(row, rep(NA, max_length - length(row))),
						max(sapply(lg, length))
						)
)
rownames(lg_matrix) <- lg_matrix[, 1]
lg_matrix=lg_matrix[,-1]
### 读入GO功能富集结果文件 
lnci_lncj_func <- fread("GBM_function_fdr_01_01_1.txt", header=F,data.table=F)
lnci_lncj_func<-as.matrix(lnci_lncj_func)
################################################
### 网络 ########################################
library(igraph)
### 读入蛋白质网络 
mcp <- fread("J:/xieyunjin/test-GBM-model/data/maxcomponent7.28.txt", header=F,data.table=F)
PPI <- graph.data.frame(mcp, directed=F)
### 读入随机网络,并转换成igraph图对象
random_net_names <- dir("J:/xieyunjin/test-GBM-model/data/randomnetwork")
random_nets <- lapply(random_net_names, function(x){
												path <- paste("J:/xieyunjin/test-GBM-model/data/randomnetwork", x, sep="/")
												rn <- fread(path, header = F,data.table=F)
												graph.data.frame(rn, directed = F)
										}
					)
### 网络顶点名
vertice_names <- get.data.frame(PPI, what = "vertices")[,1]
################################################
### 主程序，检验模块性 ############################
NROW <- nrow(lnci_lncj_func)
### 打开输出文件
outfile <- file("GBM_fit_model_p_2.txt", "w")
for(each_row in 1 : NROW){
	row_data <- lnci_lncj_func[each_row, ]
	### row_data的第6列(模块候选基因数) 与 第3列（lncRNA对共同靶基因的个数）的比例 >= 0.01
	opp1<-as.numeric(row_data[6])
	opp2<-as.numeric(row_data[3])
	if((opp1 / opp2) >= 0.01){
		input <- as.character(row_data[c(1, 2, 4)])
		### lncRNA_i 基因集合 
		lnci <- input[1]
		lnci_genes <- lg_matrix[lnci, ]
		i_g <- lnci_genes[! is.na(lnci_genes)]
		### lncRNA_j 基因集合
		lncj <- input[2]
		lncj_genes <- lg_matrix[lncj, ]
		j_g <- lncj_genes[! is.na(lncj_genes)]
		### 候选模块 基因集合
		fc <- input[3]
		fc<-as.numeric(fc)
		fc<-as.character(fc)
		lncij_genes <- intersect(i_g, j_g)
		Gs_id <- intersect(g2g_matrix[fc, - 1 ], lncij_genes)
		#把Gs_id的symbol转成entrze id
		temp_G1<-as.matrix(Gs_id)
        colnames(temp_G1)<-"HGNC.symbol"
        temp_G2<-merge(temp_G1,gene_id,by="HGNC.symbol")
        temp_G3<-temp_G2[,2]
		Gs<-unique(temp_G3)
		Gs<-as.character(Gs)
		### 检验 Gs是否全包含在网络中
		test_Gs <- all(! is.na(match(Gs, vertice_names)))##all()检查是否全部为真，若全为真返回TURE
		if(test_Gs){
			##############################################
			### 计算 蛋白质互作网络 模块内两两基因间的距离D1是否<3
			### 以及在 蛋白质互作网络 中的cpl D2是否小于随机网络
			D1 <- shortest.paths(PPI, Gs, Gs)
			if(all(D1 < 3)){
				len_Gs <- length(Gs) 
				num_pairs <- choose(len_Gs, 2) * 2
				D2 <- sum(D1) / num_pairs 
				### 输入基因集合Gs，计算随机网络的CPL
				CPLs <- lapply(random_nets, function(r_net){
														sum(shortest.paths(r_net, Gs, Gs)) / num_pairs
														}
								       )
				### 计算 p值
				p_value <- sum(CPLs < D2) / 1000
				out_value<-c(input, p_value) 
				cat(out_value, sep = "\t", file = outfile)
					cat("\n", file = outfile)
			}
		}
	}
}
close(outfile)


















