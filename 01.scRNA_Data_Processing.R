
#Single-cell data pre-processing
library(Seurat)
library(harmony)
HGNC_PCG <- read.delim("HGNC_20240123_PCG.txt")
Data_Name <- dir("File")
Sample <- Data_Name
Data <- NULL
for (i in 1:length(Data_Name)) {
  Data[[i]] <- Read10X(data.dir = paste0("Data_RAW\\",Data_Name[i]))
  Data[[i]] <- Data[[i]][rownames(Data[[i]]) %in% HGNC_PCG$symbol,]
  Data[[i]] <- CreateSeuratObject(counts = Data[[i]],project = Sample[i],min.cells = 3,min.features = 200)
  Data[[i]][["percent.mt"]] <- PercentageFeatureSet(Data[[i]], pattern = "^MT-")
  Outlier_Cutoff <- median(Data[[i]]@meta.data$nFeature_RNA) + 3*mad(Data[[i]]@meta.data$nFeature_RNA,constant = 1)
  Data[[i]] <- subset(Data[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < Outlier_Cutoff & percent.mt < 20)
  Data[[i]] <- RenameCells(Data[[i]],new.names = paste0(Sample[i],"_",Cells(Data[[i]])))
}
Data <- merge(Data[[1]],y = Data[2:length(Data_Name)])
Data <- NormalizeData(Data) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
Data <- RunHarmony(Data, group.by.vars = "orig.ident",max.iter.harmony = 20)
Data <- FindNeighbors(Data, reduction = "harmony", dims = 1:30)
Data <- FindClusters(Data, resolution = 0.8)
Data <- RunUMAP(Data, reduction = "harmony", dims = 1:30)
saveRDS(Data,file = "Result\\Data_Clustering.RDS") 


