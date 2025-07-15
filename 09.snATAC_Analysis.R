

#1.Peak calling
library(Signac)
GSM <- dir()
for (i in 1:length(GSM)) {
  fragpath <- paste0(GSM[i],"/fragments.tsv.gz")
  total_counts <- CountFragments(fragpath)
  barcodes <- total_counts[total_counts$frequency_count > 1000, ]$CB
  frags <- CreateFragmentObject(path = fragpath, cells = barcodes)
  peaks <- CallPeaks(frags,macs2.path = ".../.conda/envs/MACS2/bin/macs2")
  peak_data <- as.data.frame(peaks)
  bed_data <- data.frame(chr = peak_data$seqnames, start = peak_data$start, end = peak_data$end)
  Chromosome <- paste0("chr",c(1:22,"X","Y"))
  bed_data <- bed_data[bed_data$chr %in% Chromosome,]
  outpath <- paste0(GSM[i],"/peaks.bed")
  write.table(bed_data, file = outpath, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}


#2.Pre-processing
library(Signac)
library(Seurat)
library(harmony)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Data <- c()
File <- dir(".../Data/RAW")
for (i in 1:length(File)) {
  Fragment <- CreateFragmentObject(path = paste0(".../Data/RAW/",File[i],"/fragments.tsv.gz"))
  Counts <- FeatureMatrix(fragments = Fragment,features = combined.peaks)
  Assay <- CreateChromatinAssay(Counts, sep = c("-", "-"), fragments = Fragment, min.cells = 10,min.features = 200)
  Data[[i]] <- CreateSeuratObject(Assay, assay = "Peaks")
  Data[[i]]$Sample <- File[i]
  Annotation(Data[[i]]) <- annotations
  Data[[i]] <- NucleosomeSignal(object = Data[[i]])
  Data[[i]] <- TSSEnrichment(object = Data[[i]], fast = F)
  Data[[i]]$blacklist_ratio <- FractionCountsInRegion(Data[[i]],assay = 'Peaks',regions = blacklist_hg38_unified)
  Data[[i]] <- subset(x = Data[[i]],subset = nucleosome_signal < 4 & TSS.enrichment > 1 &nCount_Peaks > 3000 & nCount_Peaks < 30000 & blacklist_ratio < 0.01)
  Gene_Activity <- GeneActivity(Data[[i]])
  Data[[i]][['RNA']] <- CreateAssayObject(counts = Gene_Activity)
  Data[[i]] <- NormalizeData(Data[[i]],assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = median(Data[[i]]$nCount_RNA))
}
Combined <- merge(x = Data[[1]], y = Data[2:length(Data)],add.cell.ids = File)
Combined <- RunTFIDF(Combined)
Combined <- FindTopFeatures(Combined, min.cutoff = 'q0')
Combined <- RunSVD(Combined)
Combined <- RunHarmony(Combined,group.by.vars = 'Sample',reduction.use = 'lsi',assay.use = 'Peaks',project.dim = F)
Combined <- FindNeighbors(Combined, reduction = 'harmony', dims = 2:30)
Combined <- FindClusters(Combined, verbose = FALSE, algorithm = 3)
Combined <- RunUMAP(Combined, dims = 2:30, reduction = 'harmony')
saveRDS(Combined,".../Data_Processing/Data.RDS")


#3.Annotation
library(Signac)
library(Seurat)
Data <- readRDS(".../Data_Processing/Data.RDS")
SC <- readRDS(".../Reference/Reference.RDS")
Anchors <- FindTransferAnchors(reference = SC,query = Data,reduction = 'cca')
Predicted_Labels <- TransferData(anchorset = Anchors,refdata = SC$Celltype,weight.reduction = Data[['lsi']],dims = 2:30)
Data <- AddMetaData(object = Data, metadata = Predicted_Labels)
for(i in levels(Data)) {
  cells_to_reid <- WhichCells(Data, idents = i)
  newid <- names(which.max(table(Data$predicted.id[cells_to_reid])))
  Idents(Data, cells = cells_to_reid) <- newid
}
Data$Celltype <- Idents(Data)
saveRDS(Data,".../Annotation/Data.RDS")


#4.Fibroblast Integration
library(Signac)
library(Seurat)
library(ggplot2)
library(harmony)
GSE214085 <- readRDS(".../Annotation/GSE214085.RDS")
GSE195460 <- readRDS(".../Annotation/GSE195460.RDS")
GSE151302 <- readRDS(".../Annotation/GSE151302.RDS")
GSE270788 <- readRDS(".../Annotation/GSE270788.RDS")
GSE214085 <- subset(GSE214085,idents = "Fibroblast")
GSE195460 <- subset(GSE195460,idents = "Fibroblast")
GSE151302 <- subset(GSE151302,idents = "Fibroblast")
GSE270788 <- subset(GSE270788,idents = "Fibroblast")
GSE214085$Dataset <- "GSE214085"
GSE195460$Dataset <- "GSE195460"
GSE151302$Dataset <- "GSE151302"
GSE270788$Dataset <- "GSE270788"
Combined <- merge(x = GSE214085,y = list(GSE195460,GSE151302,GSE270788))
DefaultAssay(Combined) <- 'Peaks'
Combined <- RunTFIDF(Combined)
Combined <- FindTopFeatures(Combined, min.cutoff = 10)
Combined <- RunSVD(Combined)
Combined <- RunHarmony(Combined,group.by.vars = "Dataset",max_iter = 20,reduction.use = 'lsi',assay.use = 'Peaks',project.dim = F)
Combined <- FindNeighbors(Combined, reduction = 'harmony', dims = 2:30)
Combined <- FindClusters(Combined, verbose = FALSE, algorithm = 3)
Combined <- RunUMAP(Combined, dims = 2:30, reduction = 'harmony')
saveRDS(Combined,".../Combined_Fib/Combined_Fib.RDS")


#5.Fibroblast Annotation
Fib_ATAC <- readRDS(".../Combined_Fib/Combined_Fib.RDS")
Fib_SC <- readRDS(".../Reference/Fib_Reference.RDS")
Anchors <- FindTransferAnchors(reference = Fib_SC,query = Fib_ATAC,reduction = 'cca')
Predicted_Labels <- TransferData(anchorset = Anchors,refdata = Fib_SC$Subtype,weight.reduction = Fib_ATAC[['lsi']],dims = 2:30)
Fib_ATAC <- AddMetaData(object = Fib_ATAC, metadata = Predicted_Labels)
for(i in levels(Fib_ATAC)) {
  cells_to_reid <- WhichCells(Fib_ATAC, idents = i)
  newid <- names(which.max(table(Fib_ATAC$predicted.id[cells_to_reid])))
  Idents(Fib_ATAC, cells = cells_to_reid) <- newid
}
Fib_ATAC$Subtype <- Idents(Fib_ATAC)
saveRDS(Fib_ATAC,".../Combined_Fib/Fib_ATAC_Annotation.RDS")


#6.DAR Analysis
library(Signac)
library(Seurat)
library(presto)
Fib_ATAC <- readRDS(".../Combined_Fib/Fib_ATAC_Annotation.RDS")
DefaultAssay(Fib_ATAC) <- 'Peaks'
DA_Peaks <- FindAllMarkers(Fib_ATAC,min.pct = 0.1,only.pos = TRUE)
Marker <- DA_Peaks$gene
Closest_Genes <- ClosestFeature(Fib_ATAC, regions = Marker)


