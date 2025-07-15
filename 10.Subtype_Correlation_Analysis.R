

library(dplyr)
library(RColorBrewer)
library(pheatmap)
Meta_All <- read.delim(".../Meta_All.txt")
Meta_All$Subtype <- paste0(Meta_All$Celltype,"_",Meta_All$Subtype)
Meta_Endothelial <- Meta_All[Meta_All$Celltype == "Endothelial",c("orig.ident","Subtype")]
Meta_Macrophage <- Meta_All[Meta_All$Celltype == "Macrophage",c("orig.ident","Subtype")]
Meta_Tcell <- Meta_All[Meta_All$Celltype == "Tcell",c("orig.ident","Subtype")]
Meta_Bcell <- Meta_All[Meta_All$Celltype == "Bcell",c("orig.ident","Subtype")]
Meta_Fib <- Meta_All[Meta_All$Celltype == "Fib",c("orig.ident","Subtype")]
Freq_Endothelial <- data.frame(table(Meta_Endothelial)/rowSums(table(Meta_Endothelial)))
Freq_Macrophage <- data.frame(table(Meta_Macrophage)/rowSums(table(Meta_Macrophage)))
Freq_Tcell <- data.frame(table(Meta_Tcell)/rowSums(table(Meta_Tcell)))
Freq_Bcell <- data.frame(table(Meta_Bcell)/rowSums(table(Meta_Bcell)))
Freq_Fib <- data.frame(table(Meta_Fib)/rowSums(table(Meta_Fib)))
Freq_All <- rbind(Freq_Endothelial,Freq_Macrophage,Freq_Tcell,Freq_Bcell,Freq_Fib)
Result_Freq <- matrix(NA,nrow = length(unique(Meta_All$orig.ident)),ncol = length(unique(Meta_All$Subtype)),dimnames = list(unique(Meta_All$orig.ident),unique(Meta_All$Subtype)))
for (i in 1:ncol(Result_Freq)) {
  Subtype_Tem <- colnames(Result_Freq)[i]
  Freq_Tem <- Freq_All[Freq_All$Subtype == Subtype_Tem,]
  Result_Freq[match(Freq_Tem$orig.ident,rownames(Result_Freq)),i] <- Freq_Tem$Freq
}
Result_Cor <- c()
for (i in 1:ncol(Result_Freq)) {
  Tem1 <- colnames(Result_Freq)[i]
  for (j in 1:ncol(Result_Freq)) {
    Tem2 <- colnames(Result_Freq)[j]
    Cor_Data <- data.frame(Tem1 = Result_Freq[,i],Tem2 = Result_Freq[,j])
    Cor_Data <- Cor_Data[!is.na(Cor_Data[,1]),]
    Cor_Data <- Cor_Data[!is.na(Cor_Data[,2]),]
    Cor <- cor.test(Cor_Data[,1],Cor_Data[,2])$estimate
    Result_Cor <- rbind(Result_Cor,c(Tem1,Tem2,Cor))
  }
}
colnames(Result_Cor) <- c("Subtype1","Subtype2","Cor")
Result_Matrix <- matrix(as.numeric(Result_Cor[,3]),ncol = 31)
rownames(Result_Matrix) <- colnames(Result_Freq)
colnames(Result_Matrix) <- colnames(Result_Freq)
pheatmap(Result_Matrix,show_colnames = T,
         cutree_cols = 6,cutree_rows = 6,
         fontsize_row = 10,fontsize_col = 10,
         cluster_cols = TRUE,cluster_rows = TRUE,
         treeheight_row = 60,treeheight_col = 60,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")



