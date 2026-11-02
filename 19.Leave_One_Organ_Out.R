Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

dir.create(file.path(Repo_Dir, "Results/Leave_One_Organ_Out/Single_Organ"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(Repo_Dir, "Results/Leave_One_Organ_Out"))

rm(list = setdiff(ls(), "Repo_Dir"))
library(Seurat)
library(tidyverse)
library(harmony)
Lung <- readRDS(file.path(Repo_Dir, "data/organs/Lung_JoinLayers_Celltype.RDS"))
Liver <- readRDS(file.path(Repo_Dir, "data/organs/Liver_JoinLayers_Celltype.RDS"))
Kidney <- readRDS(file.path(Repo_Dir, "data/organs/Kidney_JoinLayers_Celltype.RDS"))
Heart <- readRDS(file.path(Repo_Dir, "data/organs/Heart_JoinLayers_Celltype.RDS"))
Skin <- readRDS(file.path(Repo_Dir, "data/organs/Skin_JoinLayers_Celltype.RDS"))
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
Kidney_Fib@assays[["RNA"]]@layers[["scale.data.1"]] <- NULL
Liver_Fib@assays[["RNA"]]@layers[["scale.data.1"]] <- NULL
rm(Lung, Liver, Kidney, Heart, Skin)
gc()
Fib_List <- list(Lung = Lung_Fib, Liver = Liver_Fib, Kidney = Kidney_Fib, Heart = Heart_Fib, Skin = Skin_Fib)
Gene_Del <- read.delim(file.path(Repo_Dir, "data/gene_sets/Gene_Del.txt"))

rm(Lung_Fib, Liver_Fib, Kidney_Fib, Heart_Fib, Skin_Fib)
gc()
save(Fib_List,Gene_Del,file = "Leave_Out.RData")

rm(list = setdiff(ls(), "Repo_Dir"))
library(Seurat)
library(tidyverse)
library(harmony)
load("Leave_Out.RData")
for (Leave_Out in names(Fib_List)) {
  message("Held-out organ: ", Leave_Out)
  set.seed(123)
  Fib_Train <- Fib_List[names(Fib_List) != Leave_Out]
  Fibroblast <- merge(Fib_Train[[1]], y = Fib_Train[-1])
  Fibroblast <- JoinLayers(Fibroblast)
  Fibroblast <- Fibroblast[!rownames(Fibroblast) %in% Gene_Del$Gene_Del, ]

  print(table(Fibroblast$Tissue))
  Fibroblast <- NormalizeData(Fibroblast) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  Fibroblast <- RunHarmony(Fibroblast, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
  Fibroblast <- FindNeighbors(Fibroblast, reduction = "harmony", dims = 1:30)
  Fibroblast <- FindClusters(Fibroblast, resolution = 0.8)
  Fibroblast <- RunUMAP(Fibroblast, reduction = "harmony", dims = 1:30)
  Fibroblast@misc$Leave_Out <- Leave_Out
  Fibroblast@misc$Training_Tissues <- names(Fib_Train)
  saveRDS(Fibroblast, paste0("Fibroblast_Without_", Leave_Out, "_Gene_Del.RDS"))
  rm(Fibroblast, Fib_Train)
  gc()
}

rm(list = setdiff(ls(), "Repo_Dir"))
library(Seurat)
library(tidyverse)
library(harmony)
load("Leave_Out.RData")
for (Tissue in names(Fib_List)) {
  message("Organ: ", Tissue)
  set.seed(123)
  Fibroblast <- Fib_List[[Tissue]]
  DefaultAssay(Fibroblast) <- "RNA"
  Fibroblast <- JoinLayers(Fibroblast, assay = "RNA")
  Fibroblast <- Fibroblast[!rownames(Fibroblast) %in% Gene_Del$Gene_Del, ]
  Fibroblast <- NormalizeData(Fibroblast) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  Fibroblast <- RunHarmony(Fibroblast, group.by.vars = c("orig.ident"), max.iter.harmony = 20)
  Fibroblast <- FindNeighbors(Fibroblast, reduction = "harmony", dims = 1:30)
  Fibroblast <- FindClusters(Fibroblast, resolution = 0.8)
  Fibroblast <- RunUMAP(Fibroblast, reduction = "harmony", dims = 1:30)
  saveRDS(Fibroblast, paste0("Single_Organ/Fibroblast_", Tissue, "_Gene_Del.RDS"))
}
