


library(Seurat)
library(dplyr)
library(spacexr)

#1.Reference
SC <- readRDS(".../SC.RDS")
Count_SC <- as.matrix(GetAssayData(SC, assay = "RNA", layer = "counts"))
CellType_SC <- Idents(SC)
Object_SC <- Reference(Count_SC, CellType_SC, n_max_cells = 100)

#2.ST Deconvolution
Sample_Info <- read.delim(".../ST_Sample_Info.txt")
HGNC_PCG <- read.delim(".../HGNC_20240123_PCG.txt")
Gene_Del <- read.delim(".../Blacklist.txt")
for (i in 1:nrow(Sample_Info)) {
  Dir <- paste0("RAW/",Sample_Info$Name[i])
  ST_Tem <- Load10X_Spatial(data.dir = Dir,filename = "filtered_feature_bc_matrix.h5",assay = "Spatial")
  colnames(ST_Tem) <- paste0(Sample_Info$Name[i],"_",colnames(ST_Tem))
  ST_Tem$Sample <- Sample_Info$Name[i]
  ST_Tem <- ST_Tem[rownames(ST_Tem) %in% HGNC_PCG$symbol,]
  ST_Tem <- ST_Tem[!rownames(ST_Tem) %in% Gene_Del$Gene_Del, ]
  ST_Tem <- subset(ST_Tem,subset = nFeature_Spatial > 300 & nCount_Spatial > 500)
  ST_Tem <- SCTransform(ST_Tem, assay = "Spatial", verbose = FALSE)
  Count_ST <- as.matrix(GetAssayData(ST_Tem, assay = "Spatial", layer = "counts"))
  Coords_ST <- GetTissueCoordinates(ST_Tem,scale = NULL)[,-3]
  Object_ST <- SpatialRNA(Coords_ST, Count_ST)
  RCTD <- create.RCTD(Object_ST, Object_SC,max_cores = 8)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  ST_Tem[["RCTD"]] <- CreateAssayObject(t(RCTD@results$weights))
  Path_Result <- paste0("...\\RCTD\\ST\\",Sample_Info$Name[i],"_ST.RDS")
  saveRDS(ST_Tem,file = Path_Result)
}



