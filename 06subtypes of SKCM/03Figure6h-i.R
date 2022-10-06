#=====================加个热图，两类  6 IC-lncRNAs的表达 加上age, gender or cancer stage======================
setwd("J:/xieyunjin/new_result/ConsensusClusterPlus/SKCM_new")

lnc <- read.table("SKCM_filter_lnc.txt",sep='\t',header = F,as.is = T)
lnc <- as.matrix(lnc)

symbol <- read.table("J:/xieyunjin/lncRNA_synergisty_project/pan_cancer_symbol/symbol.txt",sep='\t',header = T,as.is = T)
head(symbol)
lnc_symbol <- symbol[match(lnc,symbol$ensg),]$symbol

cluster <- read.table("cluster2.txt",sep='\t',header = F,as.is = T)
cluster$V2 <- factor(cluster$V2,levels = c(1,2))
cluster2 <- cluster[order(cluster$V2),]
cluster2$V1 <- apply(as.matrix(cluster2$V1),1,function(x) paste(unlist(strsplit(x,"\\."))[3:5],collapse = "-"))
colnames(cluster2) <- c("sample","group")

exp <- read.table("J:/xieyunjin/lncRNA_expression/SKCM_regression_1_lncRNA.txt",sep = "\t",header = T,as.is = T)
exp2<-as.matrix(exp[,-1])
rownames(exp2) <- substr(as.matrix(exp[,1]),1,15)
colnames(exp2) <- apply(as.matrix(colnames(exp2)),1,function(x) paste(unlist(strsplit(x,"\\."))[3:5],collapse = "-"))
lnc_exp <- exp2[lnc,cluster2$sample]

#age
cli0<-read.table("J:/xieyunjin/UCSC_clinic_phenotype/TCGA-SKCM.GDC_phenotype.tsv",sep='\t',header=T,as.is = T,fill=T,quote = "")
cli <- unique(cli0[cli0$sample_type.samples != "Solid Tissue Normal",c(7,8,16,21,22,29:31,49,61,65,66,85,89:91)])
write.table(cli,"J:/xieyunjin/UCSC_clinic_phenotype/SKCM_phenotype.txt",sep='\t',quote = F,row.names = F)
cli2<-merge(cli,cluster2,by.x = "submitter_id",by.y='sample')
#=====年龄差距=====
pvalue<-wilcox.test(cli2[cli2$group==1,"age_at_index.demographic"],
                    cli2[cli2$group==2,"age_at_index.demographic"],alternative = "two.sided")$p.value
res<-cbind("age","wilcox",pvalue)
write.table(res,"SKCM_subtype_clinical_test.txt",sep='\t',append = T,row.names = F,col.names = F,quote = F)
#====性别差距======
g1<-table(cli2[cli2$group==1,"gender.demographic"])
g2<-table(cli2[cli2$group==2,"gender.demographic"])
df<-as.data.frame(rbind(g1,g2))
pvalue<-fisher.test(df,alternative = "two.sided")$p.value

res2<-cbind("gender","fisher",pvalue)
write.table(res2,"SKCM_subtype_clinical_test.txt",sep='\t',append = T,row.names = F,col.names = F,quote = F)
#====pathologic_M=====
table(cli2$pathologic_M)
cli2_M <- cli2[cli2$pathologic_M !="",]
cli2_M$pathologic_M2<-cli2_M$pathologic_M
cli2_M[grep("M1",cli2_M$pathologic_M),]$pathologic_M2 <-"M1"

g1<-table(cli2_M[cli2_M$group==1,"pathologic_M2"])
g2<-table(cli2_M[cli2_M$group==2,"pathologic_M2"])
df<-as.data.frame(rbind(g1,g2))
pvalue<-fisher.test(df,alternative = "two.sided")$p.value

res3<-cbind("pathologic_M","fisher",pvalue)
write.table(res3,"SKCM_subtype_clinical_test.txt",sep='\t',append = T,row.names = F,col.names = F,quote = F)

df <- as.matrix(table(cli2_M$group,cli2_M$pathologic_M))
fisher.test(df)
#0.2098
#====pathologic_N=====
table(cli2$pathologic_N)
cli2_N <- cli2[cli2$pathologic_N != "",]
cli2_N$pathologic_N2<-cli2_N$pathologic_N
cli2_N$pathologic_N2[grep("N1",cli2_N$pathologic_N2)]<-"N1"
cli2_N$pathologic_N2[grep("N2",cli2_N$pathologic_N2)]<-"N2"

g1<-table(cli2_N[cli2_N$group==1,"pathologic_N2"])
g2<-table(cli2_N[cli2_N$group==2,"pathologic_N2"])
df<-as.data.frame(rbind(g1,g2))
#===NX局部区域或淋巴结无法评估======
df2<-data.frame("N0"=rowSums(df[,1,drop=F]),
                "N1_N3"=rowSums(df[,c(2:4)]))

pvalue1<-fisher.test(df2,alternative = "two.sided")$p.value

#==
df2<-data.frame("N0_N1"=rowSums(df[,c(1:2),drop=F]),
                "N2_N3"=rowSums(df[,c(3:4)]))

pvalue2<-fisher.test(df2,alternative = "two.sided")$p.value
#===
df2<-data.frame("N0_N2"=rowSums(df[,c(1:3),drop=F]),
                "N3"=rowSums(df[,4,drop=F]))

pvalue3<-fisher.test(df2,alternative = "two.sided")$p.value
res4<-rbind(cbind("N0 VS N1-N3",pvalue1),
            cbind("N0-N1 VS N2-N3",pvalue2),
            cbind("N0-N2 VS N3",pvalue3))
write.table(res4,"SKCM_subtype_clinical_test.txt",sep='\t',append = T,row.names = F,col.names = F,quote = F)
#====pathologic_T====
#TX 原发肿瘤无法评估
#Tis 原位癌
#T0 无原发肿瘤证据
table(cli2$pathologic_T)
cli2_T <- cli2[cli2$pathologic_T!= "",]
cli2_T$pathologic_T2<-cli2_T$pathologic_T
cli2_T$pathologic_T2[grep("T1",cli2_T$pathologic_T2)]<-"T1"
cli2_T$pathologic_T2[grep("T2",cli2_T$pathologic_T2)]<-"T2"
cli2_T$pathologic_T2[grep("T3",cli2_T$pathologic_T2)]<-"T3"
cli2_T$pathologic_T2[grep("T4",cli2_T$pathologic_T2)]<-"T4"

g1<-table(cli2_T[cli2_T$group==1,"pathologic_T2"])
g2<-table(cli2_T[cli2_T$group==2,"pathologic_T2"])
df<-as.data.frame(rbind(g1,g2))

df2<-data.frame("T1_T2"=rowSums(df[,c(2:3),drop=F]),
                "T3_T4"=rowSums(df[,c(4:5)]))

pvalue<-fisher.test(df2)$p.value
res<-cbind("T1-T2 VS T3-T4",pvalue)
write.table(res,"SKCM_subtype_clinical_test.txt",sep='\t',append = T,row.names = F,col.names = F,quote = F)

dataP<-data.frame(count=as.numeric(as.matrix(df2)),
                  stage=rep(colnames(df2),each=nrow(df2)),
                  group=rep(rownames(df2),time=ncol(df2)))

pdf("pathologic_T.pdf")
ggplot(dataP,aes(x=group,y=count,fill=stage))+geom_bar(stat='identity',position='fill')+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x='',y='percent')+
  theme(axis.text = element_text(size=12,color='black'))+
  theme(axis.title = element_text(size=15,color='black'))+
  theme(legend.text = element_text(size=12,color='black'))

dev.off()
#====pathologic_stage====
table(cli2$tumor_stage.diagnoses)
# 
# cli2$pathologic_stage2<-cli2$pathologic_stage
# cli2$pathologic_stage2[which(cli2$pathologic_stage2 %in% c("Stage I","Stage IA","Stage IB"))]<-"Stage I"
# cli2$pathologic_stage2[which(cli2$pathologic_stage2 %in% c("Stage II","Stage IIA","Stage IIB","Stage IIC"))]<-"Stage II"
# cli2$pathologic_stage2[which(cli2$pathologic_stage2 %in% c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC"))]<-"Stage III"
cli2_stage <- cli2[cli2$tumor_stage.diagnoses != "not reported",]

g1<-table(cli2_stage[cli2_stage$group==1,"tumor_stage.diagnoses"])
g2<-table(cli2_stage[cli2_stage$group==2,"tumor_stage.diagnoses"])
df<-as.data.frame(rbind(g1,g2))
#=======(I/II NOS Stage 0 Stage I Stage IA Stage IB)=====
df2<-data.frame(rowSums(df[,c(1:5),drop=F]),
                rowSums(df[,c(6:14),drop=F]))

pvalue1<-fisher.test(df2)$p.value
#=========(I/II NOS Stage 0 Stage I Stage IA Stage IB Stage II Stage IIA Stage IIB Stage IIC)====
df2<-data.frame(rowSums(df[,c(1:9),drop=F]),
                rowSums(df[,c(6:14),drop=F]))

pvalue2<-fisher.test(df2)$p.value
#===========(I/II NOS Stage 0 Stage I Stage IA Stage IB Stage II Stage IIA Stage IIB 
#Stage IIC Stage III  Stage IIIA Stage IIIB Stage IIIC)=====
df2<-data.frame(rowSums(df[,c(1:13),drop=F]),
                rowSums(df[,14,drop=F]))

pvalue3<-fisher.test(df2)$p.value

res<-rbind(cbind("Stage I VS Stage II Stage III Stage IV",pvalue1),
           cbind("Stage I Stage II VS Stage III Stage IV",pvalue2),
           cbind("Stage I Stage II Stage III VS Stage IV",pvalue3))

write.table(res,"SKCM_subtype_clinical_test.txt",sep='\t',append = T,row.names = F,col.names = F,quote = F)
####画复杂热图####
setwd("J:/xieyunjin/new_result/ConsensusClusterPlus/SKCM_new")
#library(pheatmap)
library(RColorBrewer)

lnc <- read.table("SKCM_filter_lnc.txt",sep='\t',header = F,as.is = T)
lnc <- as.matrix(lnc)

symbol <- read.table("J:/xieyunjin/lncRNA_synergisty_project/pan_cancer_symbol/symbol.txt",sep='\t',header = T,as.is = T)
head(symbol)
lnc_symbol <- symbol[match(lnc,symbol$ensg),]$symbol
lnc_symbol2 <- gsub("\\..*","",lnc_symbol)

cli <- read.table("J:/xieyunjin/UCSC_clinic_phenotype/SKCM_phenotype.txt",sep='\t',header = T,as.is = T)

cluster <- read.table("cluster2.txt",sep='\t',header = F,as.is = T)
cluster$V2 <- factor(cluster$V2,levels = c(1,2))
cluster2 <- cluster[order(cluster$V2),]
cluster2$V1 <- apply(as.matrix(cluster2$V1),1,function(x) paste(unlist(strsplit(x,"\\."))[3:5],collapse = "-"))
colnames(cluster2) <- c("sample","group")

cli2<-merge(cli,cluster2,by.x = "submitter_id",by.y='sample')
cli2$group <- factor(cli2$group,levels = c(1,2))
cli3 <- cli2[order(cli2$group),]
rownames(cli3) <- as.matrix(cli3$submitter_id)

exp <- read.table("J:/xieyunjin/lncRNA_expression/SKCM_regression_1_lncRNA.txt",sep = "\t",header = T,as.is = T)
exp2<-as.matrix(exp[,-1])
rownames(exp2) <- substr(as.matrix(exp[,1]),1,15)
colnames(exp2) <- apply(as.matrix(colnames(exp2)),1,function(x) paste(unlist(strsplit(x,"\\."))[3:5],collapse = "-"))

inter_sam <- intersect(cli3$submitter_id,colnames(exp2))
cli4 <- cli3[inter_sam,]
exp3 <- exp2[lnc,inter_sam]

# cli4$age_group <- "-"
# cli4[cli4$age < 56,]$age_group <- "<56"
# cli4[cli4$age >= 56,]$age_group <- ">=56"

cli4$stage <- cli4$tumor_stage.diagnoses
cli4[cli4$tumor_stage.diagnoses %in% c("stage 0"),]$stage <- "Stage 0"
cli4[cli4$tumor_stage.diagnoses %in% c("stage i","stage ia","stage ib"),]$stage <- "Stage I/IA/IB"
cli4[cli4$tumor_stage.diagnoses %in% c("stage ii","stage iia","stage iib","stage iic"),]$stage <- "Stage II/IIA/IIB/IIC"
cli4[cli4$tumor_stage.diagnoses %in% c("stage iii","stage iiia","stage iiib","stage iiic"),]$stage <- "Stage III/IIIA/IIIB/IIIC"
cli4[cli4$tumor_stage.diagnoses %in% c("stage iv"),]$stage <- "Stage IV"

write.table(cli4,"age_gender_stage_cli.txt",sep='\t',quote = F,row.names = F)

######breslow_depth_value#########
cli4 <- read.table("age_gender_stage_cli.txt",sep='\t',header = T,as.is = T,fill = T)

cli4$breslow_depth_value <- as.numeric(cli4$breslow_depth_value)
cli4_b <- cli4[!is.na(cli4$breslow_depth_value),]
pvalue <- wilcox.test(cli4_b[cli4_b$group == "1",]$breslow_depth_value,
                      cli4_b[cli4_b$group == "2",]$breslow_depth_value,conf.level = 0.95,
                      alternative = "two.sided")$p.value
#0.01172354
kruskal.test(breslow_depth_value~group,cli4_b)
#0.01167

library(ggplot2)
library(ggsignif)
library(RColorBrewer)
pdf("breslow_depth_value_boxplot.pdf")

cli4_b$group <- as.character(cli4_b$group)

ggplot(cli4_b,aes(x = group, y = breslow_depth_value,fill = group))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA,width = 0.5)+
  geom_jitter(aes(color= group),width = 0.2,size=1.5)+
  theme_classic()+
  theme(text=element_text(size=15,color="black"),
        legend.position = "none")+
  labs(x="",y="Breslow depth")+
  geom_signif(comparisons = list(c("1","2")),
              test = "wilcox.test",y_position = 35,
              map_signif_level =  c("***"=0.001, "**"=0.01, "*"=0.05))+
  scale_fill_manual(values = c('1'=brewer.pal(9,"Set1")[3],'2'=brewer.pal(9,"Set1")[5]))+
  scale_color_manual(values = c('1'=brewer.pal(9,"Set1")[3],'2'=brewer.pal(9,"Set1")[5]))+
  theme(plot.margin = unit(c(1,3,1,3),"cm"))

dev.off()

#######malignant_neoplasm_mitotic_count_rate############
table(cli4$malignant_neoplasm_mitotic_count_rate)
mitotic <- cli4[!is.na(cli4$malignant_neoplasm_mitotic_count_rate),]

table(mitotic$malignant_neoplasm_mitotic_count_rate)
mitotic$malignant_neoplasm_mitotic_count_rate <- as.numeric(mitotic$malignant_neoplasm_mitotic_count_rate)
wilcox.test(mitotic[mitotic$group == "1",]$malignant_neoplasm_mitotic_count_rate,
            mitotic[mitotic$group == "2",]$malignant_neoplasm_mitotic_count_rate,conf.level = 0.95,
            alternative = "two.sided")$p.value
#0.02170888
library(ggplot2)
library(ggsignif)
mitotic$group <-as.character(mitotic$group)

pdf("Mitotic_rate_boxplot.pdf")
ggplot(mitotic,aes(x=factor(group),y=malignant_neoplasm_mitotic_count_rate,fill=group))+
  geom_boxplot(alpha = 0.8,outlier.shape = NA,width = 0.5)+
  geom_jitter(aes(color=group),width = 0.2,size = 1.5)+
  theme_classic()+
  theme(text=element_text(size=15,color="black"),
        legend.position = "none")+
  labs(x="",y="Mitotic rate")+
  geom_signif(comparisons = list(c("1","2")),
              test = "wilcox.test",y_position = 30,
              map_signif_level =  c("***"=0.001, "**"=0.01, "*"=0.05))+
  scale_fill_manual(values = c('1'=brewer.pal(9,"Set1")[3],'2'=brewer.pal(9,"Set1")[5]))+
  scale_color_manual(values = c('1'=brewer.pal(9,"Set1")[3],'2'=brewer.pal(9,"Set1")[5]))+
  theme(plot.margin = unit(c(1,3,1,3),"cm"))
dev.off()

# age <- matrix(cli5$age,nrow=1)
# colnames(age) <- cli5$sample
# annotation_col = data.frame(
#   Group = factor(cli5$group))
# rownames(annotation_col) <- cli5$sample
# 
# ann_colors_age = list(
#   Group = c('1'=brewer.pal(9,"Set1")[3],'2'=brewer.pal(9,"Set1")[5]))
# 
# pdf("age_heatmap.pdf")
# pheatmap(age,cluster_rows = F,cluster_cols = F,annotation_col =annotation_col,
#          annotation_colors = ann_colors_age,show_rownames=F,show_colnames=F)
# dev.off()
# 
# annotation_col = data.frame(
#   Group = factor(cli5$group),
#   Gender = factor(cli5$gender),
#   Stage = factor(cli5$stage)
# )
# rownames(annotation_col) = cli5$sample
# 
# ann_colors = list(
#   Group = c('1'=brewer.pal(9,"Set1")[3],'2'=brewer.pal(9,"Set1")[5]),
#   Gender = c("FEMALE" = "#76B7B2", "MALE"= "#EDC948"),
#   Stage = c("Stage 0" = "#A0CBE8", "Stage I" = "#FFBE7D","Stage II" = "#8CD17D",
#             "Stage III" = "#B6992D","Stage IV" = "#86BCB6")
# )
# 
# pdf("complex_heatmap_exp_clinic.pdf")
# pheatmap(exp4,cluster_rows = F,cluster_cols = F,show_colnames=F,annotation_col =annotation_col,
#          annotation_colors = ann_colors,color = c(colorRampPalette(c("#F8F8F8","#FF7856"))(80),
#                                                   colorRampPalette(c("#FF7856","#F5725E","#BB173A","#AE123A"))(20)),
#          labels_row = lnc_symbol,cellheight=35)
# dev.off()

########试一下library(ComplexHeatmap)，不用拼图########
library(ComplexHeatmap)
#连续变量颜色
library(circlize)
library(RColorBrewer)
cli4[cli4$gender.demographic == "female",]$gender.demographic <- "Female"
cli4[cli4$gender.demographic == "male",]$gender.demographic <- "Male"

col_fun2 <- colorRamp2(
  c(18, 58, 87),  #根据值的范围设置
  c("#ff7f00", "white", "#1f78b4")
)

ha <- HeatmapAnnotation(
  Group = cli4$group,
  Age = cli4$age_at_index.demographic,
  Gender = cli4$gender.demographic,
  Stage = cli4$stage,
  col = list( 
    Age = col_fun2 , #连续
    Group = c('1'=brewer.pal(9,"Set1")[3],'2'=brewer.pal(9,"Set1")[5]),#分类
    Gender = c("Female" = "#76B7B2", "Male"= "#EDC948"),
    Stage = c("i/ii nos" = "#E6E6E6","not reported" = "#E6E6E6",
              "Stage 0" = "#A0CBE8", "Stage I/IA/IB" = "#FFBE7D","Stage II/IIA/IIB/IIC" = "#8CD17D",
              "Stage III/IIIA/IIIB/IIIC" = "#B6992D","Stage IV" = "#86BCB6"
              )
  )
)

pdf("complex_heatmap_exp_clinic_cluster2.pdf",width = 10)
Heatmap(
  exp3, 
  top_annotation = ha,row_labels =  lnc_symbol2,show_column_names = F,cluster_rows = T,cluster_columns = T,
  heatmap_legend_param = list(title = "Expression"),
  col = c(colorRampPalette(c("white","#FF7856"))(80),
            colorRampPalette(c("#FF7856","#F5725E","#BB173A","#AE123A"))(20))
)
dev.off()

pdf("complex_heatmap_exp_clinic_no_cluster2.pdf",width = 10)
Heatmap(
  exp3, 
  top_annotation = ha,row_labels =  lnc_symbol2,show_column_names = F,cluster_rows = F,cluster_columns = F,
  heatmap_legend_param = list(title = "Expression"),
  col = c(colorRampPalette(c("white","#FF7856"))(80),
          colorRampPalette(c("#FF7856","#F5725E","#BB173A","#AE123A"))(20))
)
dev.off()