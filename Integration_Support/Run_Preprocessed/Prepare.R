Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({library(Seurat); library(Matrix); library(dplyr)})
setwd(file.path(Repo_Dir, "Integration_Support/Run_Preprocessed"))
future::plan("sequential")
set.seed(123)
Fibroblast <- readRDS(file.path(Repo_Dir, "data/integration/Fibroblast_Preprocessed.RDS"))
DefaultAssay(Fibroblast) <- "RNA"
stopifnot(ncol(Fibroblast) == 507934, length(unique(Fibroblast$orig.ident)) == 680, min(table(Fibroblast$orig.ident)) >= 50, !anyDuplicated(colnames(Fibroblast)))
stopifnot(all(c("counts", "data") %in% Layers(Fibroblast[["RNA"]])), length(VariableFeatures(Fibroblast)) == 2000, ncol(Embeddings(Fibroblast, "pca")) >= 30)
Sample_Info <- read.delim(file.path(Repo_Dir, "data/metadata/sample_info.tsv"), fileEncoding = "GB18030", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(Sample_Info$Sample_ID))
Batch_Info <- Fibroblast[[]] %>% mutate(Dataset_Type = Sample_Info$Dataset_Type[match(orig.ident, Sample_Info$Sample_ID)]) %>% count(orig.ident, Tissue, Dataset_Type, Fibrosis_State, name = "Cells")
stopifnot(!anyNA(Batch_Info$Dataset_Type), !anyDuplicated(Batch_Info$orig.ident))
Reference_Info <- Batch_Info %>% arrange(desc(Cells), orig.ident) %>% group_by(Tissue, Dataset_Type, Fibrosis_State) %>% slice_head(n = 1) %>% ungroup() %>% arrange(desc(Cells)) %>% slice_head(n = 20)
write.table(Batch_Info, "Tables/Batch_Info.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Reference_Info, "Tables/Reference_Batches.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
Features <- VariableFeatures(Fibroblast)
writeLines(Features, "Tables/Integration_Features.txt")
saveRDS(list(cells = colnames(Fibroblast), genes = rownames(Fibroblast), features = Features, counts_nnz = length(LayerData(Fibroblast, layer = "counts")@x), data_nnz = length(LayerData(Fibroblast, layer = "data")@x)), "Checkpoints/Input_Info.rds", compress = FALSE)
Fibroblast[["RNA"]]$scale.data <- NULL
Thin <- CreateSeuratObject(counts = LayerData(Fibroblast, layer = "counts")[Features, , drop = FALSE], meta.data = Fibroblast[[]])
Thin[["RNA"]]$data <- LayerData(Fibroblast, layer = "data")[Features, , drop = FALSE]
VariableFeatures(Thin) <- Features
Thin[["pca"]] <- Fibroblast[["pca"]]
saveRDS(Thin, "Checkpoints/Common.rds.tmp", compress = FALSE)
stopifnot(file.rename("Checkpoints/Common.rds.tmp", "Checkpoints/Common.rds"))
cat("Prepared: existing data/HVG/PCA reused; cells", ncol(Thin), "; batches", nrow(Batch_Info), "; references", nrow(Reference_Info), "\n")
writeLines(capture.output(sessionInfo()), "Logs/prepare_sessionInfo.txt")
