Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

log_stage("load original preprocessed RNA for this method")
rm(Fibroblast); gc()
Result <- readRDS(paste0("Checkpoints/", Stage, "_Result.rds"))
Fibroblast <- readRDS(file.path(Repo_Dir, "data/integration/Fibroblast_Preprocessed.RDS"))
Input_Info <- readRDS("Checkpoints/Input_Info.rds")
stopifnot(identical(colnames(Fibroblast), Input_Info$cells), identical(rownames(Fibroblast), Input_Info$genes))
Fibroblast[["RNA"]]$scale.data <- NULL
for (Name in names(Fibroblast@graphs)) Fibroblast[[Name]] <- NULL
for (Name in names(Fibroblast@neighbors)) Fibroblast[[Name]] <- NULL
Fibroblast[[paste0("integrated.", tolower(Stage))]] <- Result$embedding
Fibroblast[[paste0("umap.", tolower(Stage))]] <- Result$umap
Fibroblast[[paste0(tolower(Stage), "_clusters")]] <- factor(Result$clusters[colnames(Fibroblast)], levels = as.character(sort(as.integer(unique(Result$clusters)))))
Idents(Fibroblast) <- paste0(tolower(Stage), "_clusters")
Fibroblast@misc$Integration <- list(method = Stage, input = "Fibroblast_Preprocessed.RDS", batch = "orig.ident", dimensions = 30, resolution = 0.8, seed = 123, retained_cells = ncol(Fibroblast), retained_samples = length(unique(Fibroblast$orig.ident)), references = if (Stage %in% c("RPCA", "CCA")) read.delim("Tables/Reference_Batches.tsv") else NULL)
Fibroblast@misc$Integration$settings <- if (Stage == "LIGER") list(algorithm = "standard iNMF", k = 30, lambda = 5, nIteration = 30, nRandomStarts = 1, nCores = 4, full_gene_library_sizes = TRUE, nonnegative_scaling = TRUE, quantile_reference = readLines("Logs/LIGER_Quantile_Reference.txt"), quantiles = 50, maxSample = 1000, seed = 123) else list(dimred = 30, knn = 20, sigma = 15, alpha = 0.1, approx = TRUE, batch_size = 5000, sketch = FALSE, seed = 123)
if (Stage == "LIGER") Fibroblast@misc$Integration$settings$scaling_implementation <- "Sparse RMS with denominator n-1; verified equivalent to rliger::scaleNotCenter"
DefaultAssay(Fibroblast) <- "RNA"
Path <- paste0("../Fibroblast_Preprocessed_", Stage, ".RDS")
log_stage(paste("save full method result", Path))
if (!file.exists(Path)) save_checkpoint(Fibroblast, Path) else log_stage("existing final file; resume reload verification")
rm(Fibroblast, Result); gc()
log_stage("reload full saved RDS and verify")
Fibroblast <- readRDS(Path)
Embedding <- Embeddings(Fibroblast, paste0("integrated.", tolower(Stage)))
UMAP <- Embeddings(Fibroblast, paste0("umap.", tolower(Stage)))
Clusters <- Fibroblast[[]][[paste0(tolower(Stage), "_clusters")]]
stopifnot(identical(colnames(Fibroblast), Input_Info$cells), identical(rownames(Fibroblast), Input_Info$genes), identical(rownames(Embedding), Input_Info$cells), identical(rownames(UMAP), Input_Info$cells), ncol(Embedding) == 30, ncol(UMAP) == 2, all(is.finite(Embedding)), all(is.finite(UMAP)), !anyNA(Clusters), identical(as.character(Idents(Fibroblast)), as.character(Clusters)))
stopifnot(identical(VariableFeatures(Fibroblast), Input_Info$features), length(LayerData(Fibroblast, layer = "counts")@x) == Input_Info$counts_nnz, length(LayerData(Fibroblast, layer = "data")@x) == Input_Info$data_nnz)
write.table(data.frame(Method = Stage, Cells = ncol(Fibroblast), Samples = length(unique(Fibroblast$orig.ident)), Clusters = nlevels(Clusters), Embedding_Dimensions = ncol(Embedding), UMAP_Dimensions = ncol(UMAP), Reload_Passed = TRUE, File = basename(Path)), paste0("Logs/", Stage, "_File_QC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
log_stage(paste("saved and reload verified; clusters", nlevels(Clusters)))
