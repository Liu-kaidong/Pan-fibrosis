

library(Seurat)
library(CytoTRACE)
Fib <- readRDS(".../Fibroblast.RDS")
Cytotrace_Result <- CytoTRACE(as.matrix(Fib[["RNA"]]$counts), ncores = 1, enableFast = T)



