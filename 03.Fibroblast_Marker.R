Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)
library(Seurat)
library(ggplot2)
library(presto)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
Fib <- readRDS(file.path(Repo_Dir, "data/fibroblast/Fibroblast_Gene_Del_Cluster_Del.RDS"))
Fib_Marker <- FindAllMarkers(Fib,only.pos = TRUE,min.pct = 0.25)
Fib_Marker <- Fib_Marker %>% filter(p_val_adj < 0.05) # Seurat Bonferroni-adjusted P values
Top50 <- Fib_Marker %>%  group_by(cluster) %>%  top_n(n = 50, wt = avg_log2FC)
#IntraDEG
Idents(Fib) <- "Subtype1"
Subtype_DEG <- list()
for (Subtype in levels(Idents(Fib))) {
  Subtype_DEG[[Subtype]] <- FindMarkers(Fib, ident.1 = "Fibrosis", ident.2 = "Normal", group.by = "Fibrosis_State", subset.ident = Subtype, min.pct = 0.25)
}
#Enrichment
GOBP <- list()
for (Subtype in unique(Top50$cluster)) {
  Marker_Genes <- Top50$gene[Top50$cluster == Subtype]
  GOBP[[as.character(Subtype)]] <- enrichGO(Marker_Genes, OrgDb = "org.Hs.eg.db", ont = "BP", keyType = "SYMBOL", pAdjustMethod = "BH")
}
