

#1.Monocle2
library(Seurat)
library(monocle)
Fib <- readRDS(".../Fibroblast.RDS")
Expr_Matrix <- GetAssayData(Fib, assay = "RNA", slot = "counts")
Phenotype_Data <- Fib@meta.data
Phenotype_Data$Cell_Type <- Idents(Fib)
Feature_Data <- data.frame(gene_short_name = rownames(Fib),row.names = rownames(Fib))
Phenotype_Data <- new("AnnotatedDataFrame",data = Phenotype_Data)
Feature_Data <- new("AnnotatedDataFrame",data = Feature_Data)
CellDataSet <- newCellDataSet(Expr_Matrix, phenoData = Phenotype_Data, featureData = Feature_Data,expressionFamily = negbinomial.size())
CellDataSet <- estimateSizeFactors(CellDataSet)
ALL_Markers <- FindAllMarkers(Fib,only.pos = T,min.pct = 0.25)
Ordering_Gene <- unique(ALL_Markers$gene)
CellDataSet <- setOrderingFilter(CellDataSet,Ordering_Gene)
CellDataSet <- reduceDimension(CellDataSet, max_components = 2, method = "DDRTree")
CellDataSet <- orderCells(CellDataSet)
plot_cell_trajectory(CellDataSet, color_by = "Cell_Type")


#2.Slingshot
library(Seurat)
library(SeuratWrappers)
library(slingshot)
Fib <- readRDS(".../Fibroblast.RDS")
Reduction <- Fib@reductions$umap@cell.embeddings
Cluster <- Fib$Subtype
names(Cluster) <- colnames(Fib)
lin1 <- getLineages(Reduction, Cluster, start.clus = 'ecmFib')
crv1 <- getCurves(lin1)
Colors <- c("matFib" = "#4975A5","ecmFib" = "#BC8F8F","myoFib" = "#93A95D","apFib" = "#CD6E72")
plot(Reduction, col = Colors[match(Cluster,names(Colors))], asp = 1, pch = 16,cex = 0.1)
lines(SlingshotDataSet(lin1), lwd = 2, col = 'black',lty = 1,cex = 1.5)
lines(SlingshotDataSet(crv1), lwd = 2, col = 'black',lty = 1,cex = 1.5)


#3.Diffusion map
library(Seurat)
library(destiny)
Fib <- readRDS(".../Fibroblast.RDS")
Meta <- Fib@meta.data
harmony_data <- Embeddings(Fib, "harmony")
dm <- DiffusionMap(harmony_data)
Meta <- Meta[names(dm@eigenvec0),]
Colors <- c("matFib" = "#4975A5","ecmFib" = "#BC8F8F","myoFib" = "#93A95D","apFib" = "#CD6E72")
plot(dm,1:2, col = Colors[Meta$Subtype1]) + theme_void()



