Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

setwd(file.path(Repo_Dir, "Integration_Support/Run_Additional"))
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
future::plan("sequential")
set.seed(123)
Stage <- Sys.getenv("GM_STAGE", "LIGER")

log_stage <- function(x) {Line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Stage, x, sep = "\t"); cat(Line, "\n"); cat(Line, "\n", file = "Logs/Stages.tsv", append = TRUE)}
save_checkpoint <- function(x, path) {saveRDS(x, paste0(path, ".tmp"), compress = FALSE); if (file.exists(path)) unlink(path); file.rename(paste0(path, ".tmp"), path)}
log_stage("START")
Fibroblast <- readRDS("../Run_Preprocessed/Checkpoints/Common.rds")
Input_Info <- readRDS("Checkpoints/Input_Info.rds")

Features <- Input_Info$features
Reduction_Name <- paste0("integrated.", tolower(Stage))
if (file.exists(paste0("Checkpoints/", Stage, "_Embedding.rds"))) {
  Embedding <- readRDS(paste0("Checkpoints/", Stage, "_Embedding.rds"))
  log_stage("resume saved embedding")
} else {
  if (Stage == "LIGER") {
    suppressPackageStartupMessages(library(rliger))
    if (!file.exists("Checkpoints/Full_Library_Sizes.rds")) {
      log_stage("read full-gene raw counts for library sizes; no HVG-only normalization")
      Raw <- readRDS(file.path(Repo_Dir, "data/integration/Fibroblast_Preprocessed.RDS"))

      Library_Sizes <- Matrix::colSums(LayerData(Raw, assay = "RNA", layer = "counts"))

      save_checkpoint(Library_Sizes, "Checkpoints/Full_Library_Sizes.rds")
      rm(Raw); gc()
    }
    Library_Sizes <- readRDS("Checkpoints/Full_Library_Sizes.rds")
    Cell_Groups <- split(colnames(Fibroblast), as.character(Fibroblast$orig.ident))
    Grouped_Cells <- unlist(Cell_Groups, use.names = FALSE)
    Fibroblast <- Fibroblast[, Grouped_Cells]
    Fibroblast$orig.ident <- factor(as.character(Fibroblast$orig.ident), levels = names(Cell_Groups))
    Fibroblast@reductions <- list()
    if (file.exists("Checkpoints/LIGER_Scaled.rds")) {
      Scaled <- readRDS("Checkpoints/LIGER_Scaled.rds")
    } else {
      log_stage("non-centered sparse RMS scaling per batch, algebraically identical to scaleNotCenter; full-gene library sizes")
      Counts <- LayerData(Fibroblast, layer = "counts")[Features, , drop = FALSE]
      Scaled <- vector("list", length(Cell_Groups))
      for (i in seq_along(Cell_Groups)) {
        Cells <- Cell_Groups[[i]]
        Normalized <- Counts[, Cells, drop = FALSE]
        Normalized@x <- Normalized@x / rep.int(Library_Sizes[Cells], diff(Normalized@p))
        RMS <- sqrt(Matrix::rowSums(Normalized * Normalized) / (ncol(Normalized) - 1))
        Inverse_RMS <- numeric(length(RMS)); Inverse_RMS[RMS > 0] <- 1 / RMS[RMS > 0]
        Scaled[[i]] <- as(Matrix::Diagonal(x = Inverse_RMS) %*% Normalized, "dgCMatrix")
        dimnames(Scaled[[i]]) <- dimnames(Normalized)
        if (i %% 50 == 0) log_stage(paste("non-negative batch scaling", i, "/", length(Cell_Groups)))
      }
      Scaled <- do.call(cbind, Scaled)

      save_checkpoint(Scaled, "Checkpoints/LIGER_Scaled.rds")
      rm(Counts); gc()
    }
    Fibroblast[["RNA"]]$ligerScaleData <- Scaled
    rm(Scaled); gc()
    if (!file.exists("Checkpoints/LIGER_INMF.rds")) {
      log_stage("native runINMF: k=30 lambda=5 iterations=30 restarts=1 seed=123 cores=4")
      Fibroblast <- rliger::runINMF(Fibroblast, k = 30, lambda = 5, datasetVar = "orig.ident", layer = "ligerScaleData", reduction = "inmf", nIteration = 30, nRandomStarts = 1, seed = 123, nCores = 4, verbose = TRUE)

      save_checkpoint(Fibroblast[["inmf"]], "Checkpoints/LIGER_INMF.rds")
    } else Fibroblast[["inmf"]] <- readRDS("Checkpoints/LIGER_INMF.rds")
    H_Cells <- rownames(Embeddings(Fibroblast, "inmf"))

    Fibroblast[["inmf"]]@misc$dataset <- factor(as.character(Fibroblast$orig.ident[match(H_Cells, colnames(Fibroblast))]), levels = names(Cell_Groups))

    log_stage("native quantileNorm: batch metadata explicitly matched to embedding cell IDs; default parameters; largest batch reference")
    writeLines(names(which.max(table(Fibroblast$orig.ident))), "Logs/LIGER_Quantile_Reference.txt")
    Fibroblast <- rliger::quantileNorm(Fibroblast, reduction = "inmf", seed = 123, verbose = TRUE)
    Embedding <- Embeddings(Fibroblast, "inmfNorm")[Input_Info$cells, , drop = FALSE]
    Fibroblast[["RNA"]]$ligerScaleData <- NULL
    rm(Library_Sizes, Cell_Groups, Grouped_Cells); gc()
  }
  if (Stage == "Scanorama") {
    if (!file.exists("Checkpoints/Scanorama_Latent.tsv")) {
      log_stage("export sparse normalized HVG data, preserve original cell IDs")
      Matrix::writeMM(LayerData(Fibroblast, layer = "data")[Features, , drop = FALSE], "Checkpoints/Scanorama_Data.mtx")
      write.table(Fibroblast[[]][, "orig.ident", drop = FALSE], "Checkpoints/Scanorama_Obs.tsv", sep = "\t", quote = FALSE, col.names = NA)
      writeLines(Features, "Checkpoints/Scanorama_Genes.txt")
      log_stage("official Scanorama integrate: all cells; dims=30 knn=20 approx=True batch_size=5000 seed=123 sketch=False")
      Code <- system2("python", "Scanorama.py")

    }
    Embedding <- as.matrix(read.delim("Checkpoints/Scanorama_Latent.tsv", row.names = 1, check.names = FALSE))[Input_Info$cells, , drop = FALSE]
  }

  save_checkpoint(Embedding, paste0("Checkpoints/", Stage, "_Embedding.rds"))
  log_stage("embedding checkpoint saved")
}
if (Stage == "LIGER") {
  rm(Fibroblast); gc()
  Fibroblast <- readRDS("../Run_Preprocessed/Checkpoints/Common.rds")
}
Fibroblast[[Reduction_Name]] <- CreateDimReducObject(Embedding, key = paste0(toupper(Stage), "_"), assay = "RNA")
if (!file.exists(paste0("Checkpoints/", Stage, "_Result.rds"))) {
  if (file.exists(paste0("Checkpoints/", Stage, "_Clusters.rds"))) {
    Clusters <- readRDS(paste0("Checkpoints/", Stage, "_Clusters.rds"))
  } else {
    if (file.exists(paste0("Checkpoints/", Stage, "_SNN.rds"))) Fibroblast[["method_snn"]] <- readRDS(paste0("Checkpoints/", Stage, "_SNN.rds")) else {
      log_stage("native FindNeighbors: dims=1:30 k=20")
      Fibroblast <- FindNeighbors(Fibroblast, reduction = Reduction_Name, dims = 1:30, k.param = 20, graph.name = c("method_nn", "method_snn"))
      save_checkpoint(Fibroblast[["method_snn"]], paste0("Checkpoints/", Stage, "_SNN.rds"))
    }
    log_stage("native FindClusters: Louvain resolution=0.8 seed=123")
    Fibroblast <- FindClusters(Fibroblast, graph.name = "method_snn", resolution = 0.8, random.seed = 123, cluster.name = "Method_Cluster")
    Clusters <- setNames(as.character(Fibroblast$Method_Cluster), colnames(Fibroblast))
    save_checkpoint(Clusters, paste0("Checkpoints/", Stage, "_Clusters.rds"))
  }
  log_stage("native full-cell RunUMAP")
  Fibroblast <- RunUMAP(Fibroblast, reduction = Reduction_Name, dims = 1:30, reduction.name = "method.umap", seed.use = 123)
  Result <- list(embedding = Fibroblast[[Reduction_Name]], umap = Fibroblast[["method.umap"]], clusters = Clusters)
  save_checkpoint(Result, paste0("Checkpoints/", Stage, "_Result.rds"))
}
source("Save_Result.R", encoding = "UTF-8")
writeLines(capture.output(sessionInfo()), paste0("Logs/", Stage, "_sessionInfo.txt"))
log_stage("DONE")
