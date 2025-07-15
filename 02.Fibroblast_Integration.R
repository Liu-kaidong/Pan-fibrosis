

#1.Fibroblast Extraction
library(Seurat)
library(tidyverse)
library(harmony)
Lung <- readRDS("Lung.RDS")
Liver <- readRDS("Liver.RDS")
Kidney <- readRDS("Kidney.RDS")
Heart <- readRDS("Heart.RDS")
Skin <- readRDS("Skin.RDS")
Lung_Fib <- subset(Lung, idents = "Fibroblast")
Liver_Fib <- subset(Liver, idents = "Fibroblast")
Kidney_Fib <- subset(Kidney, idents = "Fibroblast")
Heart_Fib <- subset(Heart, idents = "Fibroblast")
Skin_Fib <- subset(Skin, idents = "Fibroblast")
Lung_Fib$Tissue <- "Lung"
Liver_Fib$Tissue <- "Liver"
Kidney_Fib$Tissue <- "Kidney"
Heart_Fib$Tissue <- "Heart"
Skin_Fib$Tissue <- "Skin"
Fibroblast <- merge(Lung_Fib, y = c(Kidney_Fib, Liver_Fib, Skin_Fib, Heart_Fib))
Fibroblast <- JoinLayers(Fibroblast)
Gene_Del <- read.delim(".../Blacklist.txt")
Fibroblast <- Fibroblast[!rownames(Fibroblast) %in% Gene_Del$Gene_Del,]
Fibroblast <- NormalizeData(Fibroblast) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
Fibroblast <- RunHarmony(Fibroblast, group.by.vars = c("orig.ident"),max.iter.harmony = 20)
Fibroblast <- FindNeighbors(Fibroblast, reduction = "harmony", dims = 1:30)
Fibroblast <- FindClusters(Fibroblast, resolution = 0.8)
Fibroblast <- RunUMAP(Fibroblast, reduction = "harmony", dims = 1:30)
saveRDS(Fibroblast,".../Fibroblast.RDS")


#2.Fibroblast Deleting
library(Seurat)
library(ggplot2)
library(dplyr)
Fib <- readRDS(".../Fibroblast.RDS") 
Markers = c('PTPRQ','WT1','NTNG1', # Podocyte
            "RYR2","TNNT2","MYH7", # Cardiomyocyte
            "ALB","TF","TTR", # Hepatocyte
            "CLDN5","PECAM1","VWF",
            "EPCAM","KRT18","KRT19",
            "NRXN1","NRXN3",'MPZ') # Neuronal cell
DotPlot(Fib, features = split(marker_genes,rep(1:6,each = 3))) + 
  RotatedAxis() + 
  theme(panel.border = element_rect(color = "black"),
        panel.spacing = unit(2,"mm"),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 60,hjust = 1),
        panel.grid = element_line(color = "grey",linetype = 3,linewidth = 0.5),
        axis.title = element_blank())
Fib <- subset(Fib,idents = c(14,16,18,19,20,22,23,25,27,29,33,34),invert = TRUE)
saveRDS(Fib,".../Fibroblast.RDS")


#3.Fibroblast Clustering
library(Seurat)
library(ggplot2)
library(dplyr)
Fib <- readRDS(".../Fibroblast.RDS") 
Fib_Aggre <- as.matrix(AggregateExpression(Fib,features = VariableFeatures(Fib))$RNA)
P <- cor(Fib_Aggre)
calc_csi <- function(pearson_cor) {
  Dim <- dim(pearson_cor)[1]
  csi <- matrix(0,nrow = Dim, ncol = Dim,dimnames = list(rownames(pearson_cor),colnames(pearson_cor)))
  for (i in seq_len(Dim)) {
    for (j in seq_len(i)) {
      Threshold <- pearson_cor[i,j] - 0.05
      P1 <- sign(pearson_cor[i,-c(i,j)] - Threshold)
      P2 <- sign(pearson_cor[j,-c(i,j)] - Threshold)
      Num <- length(which((P1 + P2) == -2))
      csi[i,j] <- Num/Dim
    }
  }
  csi <- csi + t(csi) - diag(diag(csi))
  csi
}
CSI <- calc_csi(P)
pheatmap(CSI,show_colnames = T,
         cutree_cols = 4,
         cutree_rows = 4,
         cluster_cols = TRUE,cluster_rows = TRUE,
         treeheight_row = 60,treeheight_col = 60,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")



















