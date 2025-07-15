

library(nichenetr)
library(Seurat)
library(presto)
library(tidyverse)
Fib <- readRDS(".../Fibroblast.RDS")
lr_network <- readRDS(".../NicheNet/lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS(".../NicheNet/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS(".../NicheNet/weighted_networks_nsga2r_final.rds")
lr_network <- lr_network %>% distinct(from, to)
receiver = "ecmFib"
#1.Define a set of potential ligands for the sender-agnostic approach
expressed_genes_receiver <- get_expressed_genes(receiver, Fib, pct = 0.05)
all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
#2.Define the gene set of interest
DE_table_receiver <- FindMarkers(Fib,ident.1 = "Fibrosis",ident.2 = "Normal",group.by = "Fibrosis_State",
                                 min.pct = 0.05,only.pos = T, subset.ident = "ecmFib") %>% rownames_to_column("gene")
geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
#3.Define the background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
#4.Perform NicheNet ligand activity analysis
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
best_upstream_ligands <- ligand_activities %>% top_n(10, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,"Prioritized ligands", "Ligand activity", legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank()) 

