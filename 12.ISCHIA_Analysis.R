

library(ISCHIA)
library(robustbase)
library(data.table)
library(ggplot2)
library(Seurat)
library(dplyr)
Sample_Info <- read.delim(".../ST_Sample_Info.txt")
Sample_Info$Severity_Index <- ifelse(Sample_Info$Description_B1_Mild_B2_Moderate_B3_Severe == "B1","Mild",
                                     ifelse(Sample_Info$Description_B1_Mild_B2_Moderate_B3_Severe == "B2","Moderate",
                                            ifelse(Sample_Info$Description_B1_Mild_B2_Moderate_B3_Severe == "B3","Severe","Normal")))
ST <- c()
for(i in 1:nrow(Sample_Info)){
  ST[[i]] <- readRDS(paste0(Sample_Info$Name[i],"_ST.RDS"))
  ST[[i]]$Disease_Group <- Sample_Info$Disease[i]
  ST[[i]]$Severity_Index <- Sample_Info$Severity_Index[i]
}
Combined <- merge(ST[[1]],y = ST[2:24])
names(Combined@images) <- Sample_Info$Name
RCTD <- as.matrix(t(Combined@assays[["RCTD"]]@counts))
Combined <- Combined[,colnames(Combined) %in% rownames(RCTD)]
k.values <- 1:20
wss_values <- sapply(k.values, function(k) kmeans(as.data.frame(RCTD), k, nstart = 10)$tot.withinss)
plot(k.values, wss_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K", ylab = "Total within-cluster sum of squares",
     main = "Elbow Method for Optimal K")
Combined <- Composition.cluster(Combined, RCTD, 8)
Combined_Umap <- Composition_cluster_umap(Combined, RCTD)
Embedding <- Combined_Umap[["umap.cluster.gg"]][["data"]]
Combined <- Composition.cluster(Combined, RCTD, 8)
Embedding$CompositionCluster_CC <- Combined$CompositionCluster_CC[match(rownames(Embedding), colnames(Combined))]
Idents(Combined) <- Combined$CompositionCluster_CC
ggplot(Embedding, aes(x = x, y = y)) +
  geom_point(aes(color = CompositionCluster_CC),size = 0.1,alpha = 0.8) +
  scale_color_manual(values = c("#E8797B","#B1B0E9","#D7BBB0","#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74")) +
  scale_fill_manual(values = c("#E8797B","#B1B0E9","#D7BBB0","#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +  #Í¼Àýµã´óÐ¡
  labs(x = "umap_1",y = "umap_2") +
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(0.15,"cm"),
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 15,colour = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1),"cm"),
        plot.title = element_text(size = 20,hjust = 0.5,debug = T),
        legend.text = element_text(size = 10),
        legend.key = element_blank(),
        legend.title = element_blank())



