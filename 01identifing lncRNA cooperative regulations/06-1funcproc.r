##### c5.bp.v5.1.symbols.gmt��GO.BP���ܻ��� #########
##### ����������ţ����ڴ��� #######################
##### ��������������5��С��500�Ĺ��ܼ��� ############
setwd("F://����//����1//2015new-project//lncRNA-lncRNA Эͬ���緽���ͳ���//data//Msig-GO-BP//")
gmt<-readLines("c5.bp.v6.0.symbols.gmt")#825 lines GO,BP 
g2g<-strsplit(gmt,"\t")
GLEN<-length(g2g)
for(i in 1:GLEN){
	len<-length(g2g[[i]])
	if(len>7&len<503){
		r<-c(i,g2g[[i]][c(-1,-2)])
		###### g2g: id,genes ############################
		write.table(t(r),"g2g",append=T,sep="\t",quote=F,row.names=F,col.names=F)
		###### id2function: id,function_name ###############
		write.table(t(c(i,g2g[[i]][1])),"id2function",append=T,sep="\t",quote=F,row.names=F,col.names=F)
	}
}
