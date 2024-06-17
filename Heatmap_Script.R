heatmap.frm=read.csv("Bacterial_Genus_Abundance.csv",header=T)
row.names(heatmap.frm)<-heatmap.frm$Genus
heatmap.frm<-heatmap.frm[,3:4]
heatmap_matrix<-data.matrix(heatmap.frm)
heatmap<-heatmap.2(heatmap_matrix,col= brewer.pal(9,"RdYlBu"),scale="column",margins=c(15,15),Rowv=NA,Colv=NA,trace="none",density.info="none")
