Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

log_stage("load original preprocessed RNA for this method")
rm(Fibroblast); gc()
Result <- readRDS(paste0("Checkpoints/", Stage, "_Result.rds"))
Fibroblast <- readRDS(file.path(Repo_Dir, "data/integration/Fibroblast_Preprocessed.RDS"))
Input_Info <- readRDS("Checkpoints/Input_Info.rds")

Fibroblast[["RNA"]]$scale.data <- NULL
for (Name in names(Fibroblast@graphs)) Fibroblast[[Name]] <- NULL
for (Name in names(Fibroblast@neighbors)) Fibroblast[[Name]] <- NULL
Fibroblast[[paste0("integrated.", tolower(Stage))]] <- Result$embedding
Fibroblast[[paste0("umap.", tolower(Stage))]] <- Result$umap
Fibroblast[[paste0(tolower(Stage), "_clusters")]] <- factor(Result$clusters[colnames(Fibroblast)], levels = as.character(sort(as.integer(unique(Result$clusters)))))
Idents(Fibroblast) <- paste0(tolower(Stage), "_clusters")
Fibroblast@misc$Integration <- list(method = Stage, input = "Fibroblast_Preprocessed.RDS", batch = "orig.ident", dimensions = 30, resolution = 0.8, seed = 123, retained_cells = ncol(Fibroblast), retained_samples = length(unique(Fibroblast$orig.ident)), references = if (Stage %in% c("RPCA", "CCA")) read.delim("Tables/Reference_Batches.tsv") else NULL)
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

write.table(data.frame(Method = Stage, Cells = ncol(Fibroblast), Samples = length(unique(Fibroblast$orig.ident)), Clusters = nlevels(Clusters), Embedding_Dimensions = ncol(Embedding), UMAP_Dimensions = ncol(UMAP), Reload_Passed = TRUE, File = basename(Path)), paste0("Logs/", Stage, "_File_QC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
log_stage(paste("saved and reload verified; clusters", nlevels(Clusters)))
