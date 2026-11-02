Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

setwd(file.path(Repo_Dir, "Integration_Support/Run_Preprocessed"))
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(dplyr)})
future::plan("sequential")
set.seed(123)
Stage <- Sys.getenv("GM_STAGE", "RPCA")

log_stage <- function(x) {Line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Stage, x, sep = "\t"); cat(Line, "\n"); cat(Line, "\n", file = "Logs/Stages.tsv", append = TRUE)}
save_checkpoint <- function(x, path) {saveRDS(x, paste0(path, ".tmp"), compress = FALSE); if (file.exists(path)) unlink(path); file.rename(paste0(path, ".tmp"), path)}
log_stage("START")
if (Stage %in% c("RPCA", "CCA", "scVI")) {
  Fibroblast <- readRDS("Checkpoints/Common.rds")
  Features <- readLines("Tables/Integration_Features.txt")
  Reduction_Name <- paste0("integrated.", tolower(Stage))
  if (Stage == "scVI") {Fibroblast[["RNA"]]$scale.data <- NULL; gc()}
  if (file.exists(paste0("Checkpoints/", Stage, "_Embedding.rds"))) {
    Embedding <- readRDS(paste0("Checkpoints/", Stage, "_Embedding.rds")); log_stage("resume embedding checkpoint")
  } else {
    if (Stage %in% c("RPCA", "CCA")) {
      Cell_Groups <- split(colnames(Fibroblast), as.character(Fibroblast$orig.ident))
      Reference_Info <- read.delim("Tables/Reference_Batches.tsv")
      Reference <- sort(match(Reference_Info$orig.ident, names(Cell_Groups)));
      PCA <- Fibroblast[["pca"]]
      if (file.exists(paste0("Checkpoints/", Stage, "_Anchors.rds"))) {
        Anchors <- readRDS(paste0("Checkpoints/", Stage, "_Anchors.rds"))
      } else {
        Data <- LayerData(Fibroblast, layer = "data")[Features, , drop = FALSE]
        Fib_List <- vector("list", length(Cell_Groups)); names(Fib_List) <- names(Cell_Groups)
        for (i in seq_along(Cell_Groups)) {
          Fib_List[[i]] <- CreateSeuratObject(CreateAssayObject(data = Data[, Cell_Groups[[i]], drop = FALSE]))
          VariableFeatures(Fib_List[[i]]) <- Features
          Fib_List[[i]] <- ScaleData(Fib_List[[i]], features = Features, verbose = FALSE)
          if (Stage == "RPCA") Fib_List[[i]] <- RunPCA(Fib_List[[i]], features = Features, npcs = 30, seed.use = 123, verbose = FALSE)
          if (i %% 25 == 0) log_stage(paste("sample objects", i, "/", length(Cell_Groups)))
        }
        rm(Data); gc()
        if (Stage == "CCA") {
          source("Native_Pairs.R", encoding = "UTF-8")
          Anchors <- find_anchors_pairs(Fib_List, Features, Reference, "Checkpoints/CCA_Pair_Parts", Workers = 6L, Pairs_Per_Job = 50L)
        } else {
          source("Native_Anchors.R", encoding = "UTF-8")
          Anchors <- find_anchors_chunked(Fib_List, Features, Reference, Stage, paste0("Checkpoints/", Stage, "_Anchor_Parts"), Workers = 4L, Chunk_Cells = 12000L)
        }
        save_checkpoint(Anchors, paste0("Checkpoints/", Stage, "_Anchors.rds"))
        rm(Fib_List); gc()
      }
      Anchors@reference.objects <- sort(unique(Anchors@reference.objects))
      source("Reference_Tree.R", encoding = "UTF-8")
      Sample_Tree <- reference_sample_tree(Anchors)
      save_checkpoint(Sample_Tree, paste0("Checkpoints/", Stage, "_Reference_Tree.rds"))
      log_stage("native IntegrateEmbeddings: k.weight=30; explicit default-equivalent reference tree")
      Integrated <- IntegrateEmbeddings(Anchors, reductions = PCA, new.reduction.name = Reduction_Name, dims.to.integrate = 1:30, k.weight = 30, sample.tree = Sample_Tree, verbose = TRUE)
      Embedding <- Embeddings(Integrated[[Reduction_Name]])[colnames(Fibroblast), , drop = FALSE]
      rm(Anchors, Integrated, PCA); gc()
    }
    if (Stage == "scVI") {
      if (!file.exists("Checkpoints/scVI_Latent.tsv")) {
        log_stage("export sparse raw HVG counts")
        Matrix::writeMM(LayerData(Fibroblast, layer = "counts")[Features, ], "Checkpoints/scVI_Counts.mtx")
        write.table(Fibroblast[[]][, "orig.ident", drop = FALSE], "Checkpoints/scVI_Obs.tsv", sep = "\t", quote = FALSE, col.names = NA)
        writeLines(Features, "Checkpoints/scVI_Genes.txt")
        log_stage("scVI GPU training; warm-up=50 epochs; max=200")
        Code <- system2("python", "scVI.py")

      }
      Embedding <- as.matrix(read.delim("Checkpoints/scVI_Latent.tsv", row.names = 1, check.names = FALSE))[colnames(Fibroblast), , drop = FALSE]
    }

    save_checkpoint(Embedding, paste0("Checkpoints/", Stage, "_Embedding.rds"))
    log_stage("embedding checkpoint saved")
  }
  Fibroblast[[Reduction_Name]] <- CreateDimReducObject(Embedding, key = paste0(toupper(Stage), "_"), assay = "RNA")
  if (!file.exists(paste0("Checkpoints/", Stage, "_Result.rds"))) {
    if (file.exists(paste0("Checkpoints/", Stage, "_Clusters.rds"))) {
      Clusters <- readRDS(paste0("Checkpoints/", Stage, "_Clusters.rds"))
    } else {
      if (file.exists(paste0("Checkpoints/", Stage, "_SNN.rds"))) {
        Fibroblast[["method_snn"]] <- readRDS(paste0("Checkpoints/", Stage, "_SNN.rds"))
      } else {
        log_stage("native FindNeighbors")
        Fibroblast <- FindNeighbors(Fibroblast, reduction = Reduction_Name, dims = 1:30, graph.name = c("method_nn", "method_snn"))
        save_checkpoint(Fibroblast[["method_snn"]], paste0("Checkpoints/", Stage, "_SNN.rds"))
      }
      log_stage("native FindClusters: Louvain, resolution=0.8, seed=123")
      Fibroblast <- FindClusters(Fibroblast, graph.name = "method_snn", resolution = 0.8, random.seed = 123, cluster.name = "Method_Cluster")
      Clusters <- setNames(as.character(Fibroblast$Method_Cluster), colnames(Fibroblast))
      save_checkpoint(Clusters, paste0("Checkpoints/", Stage, "_Clusters.rds"))
    }
    log_stage("native RunUMAP: all retained cells")
    Fibroblast <- RunUMAP(Fibroblast, reduction = Reduction_Name, dims = 1:30, reduction.name = "method.umap", seed.use = 123)
    Result <- list(embedding = Fibroblast[[Reduction_Name]], umap = Fibroblast[["method.umap"]], clusters = Clusters)
    save_checkpoint(Result, paste0("Checkpoints/", Stage, "_Result.rds"))
  }
}

source("Save_Result.R", encoding = "UTF-8")
writeLines(capture.output(sessionInfo()), paste0("Logs/", Stage, "_sessionInfo.txt"))
log_stage("DONE")
