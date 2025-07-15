


library(tidyverse)
library(Seurat)
library(mistyR)
future::plan(future::multisession)
Sample_Info <- read.delim(".../ST_Sample_Info.txt")
result.folders <- c()
for(i in 1:nrow(Sample_Info)){
  Sample <- Sample_Info$Name[i]
  seurat_vs <- readRDS(paste0(Sample,"_ST.RDS"))
  composition <- as.data.frame(t(seurat_vs[["RCTD"]]$data)) 
  colnames(composition) <- make.names(colnames(composition))
  Positions_Path <- paste0(".../",Sample,"/spatial/tissue_positions_list.csv")
  Position <- read.csv(Path,header = F)
  Position[,1] <- paste0(Sample,"_",Position[,1])
  geometry <- Position[match(rownames(composition), Position[,1]),c(3,4)]
  rownames(geometry) <- rownames(composition)
  colnames(geometry) <- c("row", "col")
  #MISTy
  result.folder <- create_initial_view(composition) %>% add_paraview(geometry, l = 2) %>% add_paraview(geometry, l = 5)  %>%
    run_misty(results.folder = paste0(".../MISTy/", Sample))
  names(result.folder) <- Sample
  result.folders <- c(result.folders,result.folder)
}
Sample_fibrosis <- Sample_Info$Name[Sample_Info$Disease == "Fibrosis"]
result.folders.fibrosis <- result.folders[Sample_fibrosis]
misty.results.fibrosis <- collect_results(result.folders.fibrosis)
plot_interaction_heatmap(misty.results.fibrosis, "intra", cutoff = 0.5, clean = F)
plot_interaction_heatmap(misty.results.fibrosis, "para.2", cutoff = 0.5, clean = F)
plot_interaction_heatmap(misty.results.fibrosis, "para.5", cutoff = 0.5, clean = F)



