Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

library(Seurat)
library(tidyverse)
dir.create(file.path(Repo_Dir, "Results/Integration_Annotation"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(Repo_Dir, "Results/Integration_Annotation"))
Files <- c(CCA = "2.CCA/Fibroblast_Annotation_CCA.RDS", RPCA = "3.RPCA/Fibroblast_Annotation_RPCA.RDS", scVI = "5.scVI/Fibroblast_Annotation_scVI.RDS", LIGER = "6.LIGER/Fibroblast_Annotation_LIGER.RDS", Scanorama = "7.Scanorama/Fibroblast_Annotation_Scanorama.RDS")
Order <- c("ACSM3+myoFib", "ecmFib", "ACTA2+myoFib", "apFib")
Cell_Type_Color <- c("ACSM3+myoFib" = "#4975A5", "ecmFib" = "#BC8F8F", "ACTA2+myoFib" = "#93A95D", "apFib" = "#CD6E72")
marker_genes <- c("ACSM3", "COL4A3", "MECOM", "CDH19", "THY1", "COL1A1", "COL3A1", "COL1A2", "ACTA2", "TAGLN", "MYL9", "TPM2", "CD74", "HLA-DRA", "HLA-DRB1", "SLPI")
for (Method in names(Files)) {
  Fibroblast <- readRDS(file.path(file.path(Repo_Dir, "data/integration"), Files[[Method]]))
  Fibroblast$Annotation <- gsub("\\+ +", "+", as.character(Idents(Fibroblast)))
  if (any(Fibroblast$Annotation == "myoFib")) Fibroblast$Annotation <- recode(Fibroblast$Annotation, "matFib" = "ACSM3+myoFib", "myoFib" = "ACTA2+myoFib")
  stopifnot(all(Fibroblast$Annotation %in% Order))
  Fibroblast$Annotation <- factor(Fibroblast$Annotation, levels = Order)
  Idents(Fibroblast) <- "Annotation"
  Reduction <- if (paste0("umap.", tolower(Method)) %in% Reductions(Fibroblast)) paste0("umap.", tolower(Method)) else "umap"
  stopifnot(Reduction %in% Reductions(Fibroblast))
  P1 <- DimPlot(Fibroblast, reduction = Reduction, cols = Cell_Type_Color) + labs(title = Method, color = NULL)
  P2 <- DotPlot(Fibroblast, features = marker_genes) +
    coord_flip() +
    scale_color_gradientn(colours = c("#ECEBEB", "#ECEBEB", "#DC8140", "#8C2629")) +
    labs(x = NULL, y = NULL)
  ggsave(paste0(Method, "_UMAP.pdf"), plot = P1, width = 6, height = 5)
  ggsave(paste0(Method, "_Markers.pdf"), plot = P2, width = 6, height = 7)
  write_tsv(Fibroblast[[]] %>% count(Annotation, name = "Cells") %>% mutate(Method = Method), paste0(Method, "_Cell_Counts.tsv"))
  rm(Fibroblast)
  gc()
}
