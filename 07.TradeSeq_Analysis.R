

library(Seurat)
library(SingleCellExperiment)
library(SeuratWrappers)
library(slingshot)
library(tradeSeq)
library(RColorBrewer)
Fib <- readRDS(".../Fibroblast.RDS")
sim <- as.SingleCellExperiment(Fib)
sim <- slingshot(sim, clusterLabels = 'Subtype', reducedDim = 'UMAP',  start.clus= "ecmFib", end.clus = NULL)
counts <- sim@assays@data$counts
crv <- SlingshotDataSet(sim)
set.seed(5)
icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, nGenes = 500, verbose = T) 
pseudotime <- slingPseudotime(crv, na = FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,nknots = 6, verbose = FALSE)
patternRes <- patternTest(sce)
#Visualization
patternRes$Gene <- rownames(patternRes)
patternRes <- patternRes[order(patternRes$pvalue,decreasing = F),]
patternRes <- patternRes[!is.na(patternRes$pvalue),]
patternRes$FDR <- p.adjust(patternRes$pvalue)
patternRes$log2FC <- log2(patternRes$fcMedian)
patternRes$log2FC[patternRes$log2FC > 6] <- 6
patternRes$log2FC[patternRes$log2FC < -6] <- -6
patternRes$log10FDR <- -log10(patternRes$FDR)
patternRes$curve_y <- case_when(
  patternRes$log2FC > 0 ~ 1/(patternRes$log2FC-1) + (-log10(0.05)),
  patternRes$log2FC <= 0 ~ 1/(-patternRes$log2FC-1) + (-log10(0.05)))
patternRes$group2 <- case_when(
  patternRes$`log10FDR` > patternRes$curve_y & patternRes$log2FC >= 1 ~ 'up',
  patternRes$`log10FDR` > patternRes$curve_y & patternRes$log2FC <= -1 ~ 'down',
  TRUE~ 'none')
patternRes$group2 <- factor(patternRes$group2, levels = c( "up", "down", "none"))
TF <- read.delim("Homo_sapiens_TF.txt")
Curve <- function(x){
  Pvalue <- 0.05
  LogFC <- 1
  Input <- seq(0.001, x, by = 0.001)
  Y <- 1/(Input) + (-log10(Pvalue))
  dff <- rbind(data.frame(x = Input + LogFC, y = Y),data.frame(x = -(Input + LogFC), y = Y))
  return(dff)
}
dff_curve <- Curve(10)
mycol<- c("#8A292E", "#404472", "#d8d8d8")
patternRes <- patternRes[rownames(patternRes) %in% TF$Symbol,]
patternRes <- patternRes[order(patternRes$pvalue),]
patternRes$FDR <- p.adjust(patternRes$pvalue)
patternRes$log10FDR <- -log10(patternRes$FDR)
patternRes$log10FDR[patternRes$log10FDR == "Inf"] <- 15
patternRes$Label <- patternRes$Gene
patternRes$Label[patternRes$group2 %in% "none"] <- ""
ggplot(data = patternRes,aes(x = log2FC, y = log10FDR, color = group2)) +
  geom_point(size = 6) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by = 3))+
  scale_y_continuous(expand = expansion(add = c(2, 0)),limits= c(0, 17), breaks = seq(0, 20, by = 5)) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  geom_line(data = dff_curve,aes(x = x, y = y), color= "black",lty = "dashed", size = 0.7) +
  theme_classic() +
  theme(axis.line.x = element_line(size = 1,color = "black"),
        axis.line.y = element_line(size = 1,color = "black"),
        axis.text.x = element_text(size = 15,color = "black"),
        axis.text.y = element_text(size = 15,color = "black"),
        axis.title.x = element_text(size = 18,color = "black"),
        axis.title.y = element_text(size = 18,color = "black"),
        axis.ticks = element_line(size = 1,color = "black"),axis.ticks.length = unit(0.25,"cm"),
        panel.background = element_rect(fill = "white"),
        plot.margin = unit(c(1, 1, 1, 1),"cm"),
        legend.text = element_text(size = 15),
        legend.key = element_blank(),
        legend.title = element_blank()) +
  geom_text_repel(aes(label = Label),size = 4,max.overlaps = 10000,              
                  box.padding = unit(0.5,"lines"),
                  point.padding = unit(0, "lines"))

plotSmoothers(sce, counts, "MYC", alpha = 0) + labs(title = "MYC") +
  theme_classic() + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 14),legend.text = element_text(size = 14))
plotSmoothers(sce, counts, "SMAD3", alpha = 0) + labs(title = "SMAD3") +
  theme_classic() + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 14),legend.text = element_text(size = 14))
plotSmoothers(sce, counts, "CDKN1A", alpha = 0) + labs(title = "CDKN1A") +
  theme_classic() + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 14),legend.text = element_text(size = 14))
plotSmoothers(sce, counts, "HIF1A", alpha = 0) + labs(title = "HIF1A") +
  theme_classic() + theme(axis.title = element_text(size = 15),axis.text = element_text(size = 14),legend.text = element_text(size = 14))


