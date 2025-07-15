

library(Seurat)
library(tidyverse)
library(decoupleR)
library(pheatmap)
Fib <- readRDS(".../Fibroblast.RDS")
PROGENy <- get_progeny(organism = 'human', top = 500)
Mat <- as.matrix(Fib[["RNA"]]$data)
Result <- run_mlm(mat = Mat, net = PROGENy, .source = 'source', .target = 'target',.mor = 'weight', minsize = 5)
Fib[["Pathway"]] <- Result %>% pivot_wider(id_cols = 'source', names_from = 'condition',values_from = 'score') %>% column_to_rownames('source') %>% Seurat::CreateAssayObject(.)
DefaultAssay(Fib) <- "Pathway"
Fib <- ScaleData(Fib)
Fib[["Pathway"]]$data <- Fib[["Pathway"]]$scale.data
Plot_Data <- t(as.matrix(Fib[["Pathway"]]$data)) %>% as.data.frame() %>% mutate(cluster = Idents(Fib)) %>% 
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>% group_by(cluster, source) %>% summarise(mean = mean(score))
Plot_Data1 <- Plot_Data %>% pivot_wider(id_cols = 'cluster', names_from = 'source',values_from = 'mean') %>% column_to_rownames('cluster') %>% as.matrix()
pheatmap(t(Plot_Data1), color = colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB","#4575B4")))(100),
         border_color = "black",treeheight_row = 20,treeheight_col = 20)

