Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)
find_anchors_pairs <- function(Objects, Features, Reference_Index, Cache_Dir, Workers = 4L, Pairs_Per_Job = 50L) {
  Reference_Index <- sort(unique(Reference_Index))
  Nc <- vapply(Objects, ncol, integer(1)); Offset <- c(0, head(cumsum(Nc), -1))
  Pairs <- expand.grid(seq_along(Objects), seq_along(Objects))
  Pairs <- Pairs[Pairs$Var1 < Pairs$Var2 & (Pairs$Var1 %in% Reference_Index | Pairs$Var2 %in% Reference_Index), , drop = FALSE]
  Jobs <- split(seq_len(nrow(Pairs)), ceiling(seq_len(nrow(Pairs)) / Pairs_Per_Job))
  dir.create(Cache_Dir, recursive = TRUE, showWarnings = FALSE)
  write.table(Pairs, file.path(Cache_Dir, "Pair_Order.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  Old_Max <- getOption("future.globals.maxSize"); options(future.globals.maxSize = 4 * 1024^3)
  on.exit(options(future.globals.maxSize = Old_Max), add = TRUE)
  future::plan(future::multisession, workers = Workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  Running <- list(); Done <- integer(); Next <- 1L
  say <- function(z) {if (exists("log_stage", envir = .GlobalEnv)) get("log_stage", envir = .GlobalEnv)(z) else message(z)}
  say(paste("native CCA unique reference-query pairs", nrow(Pairs), "jobs", length(Jobs)))
  while (length(Done) < length(Jobs)) {
    while (length(Running) < Workers && Next <= length(Jobs)) {
      id <- Next; Next <- Next + 1L
      File <- file.path(Cache_Dir, sprintf("Part_%04d.rds", id))
      if (file.exists(File)) {Done <- c(Done, id); next}
      Job_Pairs <- Pairs[Jobs[[id]], , drop = FALSE]
      Indices <- sort(unique(c(Job_Pairs$Var1, Job_Pairs$Var2)))
      Job_Objects <- Objects[Indices]
      Job <- future::future({
        suppressPackageStartupMessages(library(Seurat))
        future::plan(future::sequential)
        options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
        set.seed(123)
        Parts <- vector("list", nrow(job.pairs)); Prototype <- NULL
        for (j in seq_len(nrow(job.pairs))) {
          i1 <- job.pairs$Var1[j]; i2 <- job.pairs$Var2[j]
          A <- FindIntegrationAnchors(object.list = job.objects[match(c(i1, i2), job.indices)], anchor.features = job.features, scale = FALSE, reduction = "cca", dims = 1:30, k.filter = NA, reference = which(c(i1, i2) %in% job.references), verbose = TRUE)
          Forward <- A@anchors[A@anchors$dataset1 == 1 & A@anchors$dataset2 == 2, , drop = FALSE]
          Forward$dataset1 <- i1; Forward$dataset2 <- i2
          Parts[[j]] <- Forward
          if (j == 1L) {A@object.list <- lapply(A@object.list, function(x) suppressWarnings(DietSeurat(x, features = job.features[1:2]))); Prototype <- A}
        }
        list(Anchors = do.call(rbind, Parts), Prototype = Prototype)
      }, globals = list(job.objects = Job_Objects, job.indices = Indices, job.pairs = Job_Pairs, job.features = Features, job.references = Reference_Index), packages = c("Seurat", "future"), seed = 123)
      Running[[as.character(id)]] <- list(Future = Job, File = File)
      say(paste("native pair job START", id, "/", length(Jobs)))
    }
    for (id in names(Running)) {
      if (!future::resolved(Running[[id]]$Future)) next
      Part <- future::value(Running[[id]]$Future)
      saveRDS(Part, paste0(Running[[id]]$File, ".tmp"), compress = FALSE)
      file.rename(paste0(Running[[id]]$File, ".tmp"), Running[[id]]$File)
      say(paste("native pair job DONE", id, "/", length(Jobs), "forward anchors", nrow(Part$Anchors)))
      Done <- c(Done, as.integer(id)); Running[[id]] <- NULL
      rm(Part); gc()
    }
    if (length(Done) < length(Jobs)) Sys.sleep(1)
  }
  future::plan(future::sequential)
  Parts <- vector("list", length(Jobs)); First <- NULL
  for (id in seq_along(Jobs)) {Part <- readRDS(file.path(Cache_Dir, sprintf("Part_%04d.rds", id))); Parts[[id]] <- Part$Anchors; if (id == 1L) First <- Part$Prototype}
  Forward <- do.call(rbind, Parts)

  Mirrored <- Forward[, c("cell2", "cell1", "score", "dataset2", "dataset1")]; names(Mirrored) <- names(Forward)
  First@anchors <- rbind(Forward, Mirrored)
  First@offsets <- Offset; First@reference.objects <- Reference_Index
  First@object.list <- lapply(Objects, function(x) suppressWarnings(DietSeurat(x, features = Features[1:2])))
  First@command@params$reference <- Reference_Index; First@command@params$reduction <- "cca"
  methods::validObject(First)
  say(paste("native unique CCA pairs combined", nrow(First@anchors), "anchor rows"))
  First
}
