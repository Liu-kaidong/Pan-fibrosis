


library(Seurat)
library(ggplot2)
library(presto)
library(clusterProfiler)
library(org.Hs.eg.db)
Fib <- readRDS(".../Fibroblast.RDS")
Fib_Marker <- FindAllMarkers(Fib,only.pos = TRUE,min.pct = 0.25)
Top50 <- Fib_Marker %>%  group_by(subtype) %>%  top_n(n = 50, wt = avg_log2FC)
#IntraDEG
library(Seurat)
library(presto)
Fib <- readRDS(".../Fibroblast.RDS")
Subtype1_DEG <- FindMarkers(Fib,ident.1 = "Fibrosis",ident.2 = "Normal",group.by = "Fibrosis_State",subset.ident = "Subtype1",min.pct = 0.25)
Subtype2_DEG <- FindMarkers(Fib,ident.1 = "Fibrosis",ident.2 = "Normal",group.by = "Fibrosis_State",subset.ident = "Subtype2",min.pct = 0.25)
Subtype3_DEG <- FindMarkers(Fib,ident.1 = "Fibrosis",ident.2 = "Normal",group.by = "Fibrosis_State",subset.ident = "Subtype3",min.pct = 0.25)
Subtype4_DEG <- FindMarkers(Fib,ident.1 = "Fibrosis",ident.2 = "Normal",group.by = "Fibrosis_State",subset.ident = "Subtype4",min.pct = 0.25)
#Enrichment
GOBP <- enrichGO(Marker_Genes, OrgDb = "org.Hs.eg.db",ont = "BP",keyType = "SYMBOL")


