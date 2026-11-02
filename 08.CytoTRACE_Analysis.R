Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

library(Seurat)
library(CytoTRACE)
Fib <- readRDS(file.path(Repo_Dir, "data/fibroblast/Fibroblast_Gene_Del_Cluster_Del.RDS"))
Cytotrace_Result <- CytoTRACE(as.matrix(Fib[["RNA"]]$counts), ncores = 1, enableFast = T)
