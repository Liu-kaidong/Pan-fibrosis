Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)
find_anchors_chunked <- function(Objects, Features, Reference_Index, Method, Cache_Dir, Workers = 4L, Chunk_Cells = 12000L) {
  Reference_Index <- sort(unique(Reference_Index))
  dir.create(Cache_Dir, recursive = TRUE, showWarnings = FALSE)
  Nc <- vapply(Objects, ncol, integer(1))
  Offset <- c(0, head(cumsum(Nc), -1))
  Queries <- setdiff(seq_along(Objects), Reference_Index)
  Chunks <- list(); Current <- integer(); Count <- 0L
  for (i in Queries) {
    if (length(Current) && Count + Nc[i] > Chunk_Cells) {Chunks[[length(Chunks) + 1L]] <- Current; Current <- integer(); Count <- 0L}
    Current <- c(Current, i); Count <- Count + Nc[i]
  }
  if (length(Current)) Chunks[[length(Chunks) + 1L]] <- Current

  Old_Max <- getOption("future.globals.maxSize"); options(future.globals.maxSize = 4 * 1024^3)
  on.exit(options(future.globals.maxSize = Old_Max), add = TRUE)
  future::plan(future::multisession, workers = Workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  Running <- list(); Done <- integer(); Next <- 1L
  say <- function(z) {if (exists("log_stage", envir = .GlobalEnv)) get("log_stage", envir = .GlobalEnv)(z) else message(z)}
  while (length(Done) < length(Chunks)) {
    while (length(Running) < Workers && Next <= length(Chunks)) {
      id <- Next; Next <- Next + 1L
      File <- file.path(Cache_Dir, sprintf("Part_%03d.rds", id))
      if (file.exists(File)) {Done <- c(Done, id); next}
      Indices <- sort(c(Reference_Index, Chunks[[id]]))
      Job_Objects <- Objects[Indices]
      Job <- future::future({
        suppressPackageStartupMessages(library(Seurat))
        future::plan(future::sequential)
        options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
        set.seed(123)
        A <- FindIntegrationAnchors(object.list = job.objects, anchor.features = job.features, scale = FALSE, reduction = job.method, dims = 1:30, k.filter = NA, reference = job.references, verbose = TRUE)
        A@object.list <- lapply(A@object.list, function(x) suppressWarnings(DietSeurat(x, features = job.features[1:2])))
        A
      }, globals = list(job.objects = Job_Objects, job.features = Features, job.method = tolower(Method), job.references = match(Reference_Index, Indices)), packages = c("Seurat", "future"), seed = 123)
      Running[[as.character(id)]] <- list(Future = Job, Indices = Indices, File = File)
      say(paste("anchor chunk START", id, "/", length(Chunks), "query cells", sum(Nc[Chunks[[id]]])))
    }
    for (id in names(Running)) {
      if (!future::resolved(Running[[id]]$Future)) next
      A <- future::value(Running[[id]]$Future)
      Part <- list(Anchor = A, Indices = Running[[id]]$Indices)
      saveRDS(Part, paste0(Running[[id]]$File, ".tmp"), compress = FALSE)
      file.rename(paste0(Running[[id]]$File, ".tmp"), Running[[id]]$File)
      say(paste("anchor chunk DONE", id, "/", length(Chunks), "anchor rows", nrow(A@anchors)))
      Done <- c(Done, as.integer(id)); Running[[id]] <- NULL
      rm(A, Part); gc()
    }
    if (length(Done) < length(Chunks)) Sys.sleep(1)
  }
  future::plan(future::sequential)
  All <- vector("list", length(Chunks)); First <- NULL
  for (id in seq_along(Chunks)) {
    Part <- readRDS(file.path(Cache_Dir, sprintf("Part_%03d.rds", id)))
    A <- Part$Anchor@anchors
    d1 <- as.integer(A[, 4]); d2 <- as.integer(A[, 5])
    g1 <- Part$Indices[d1]; g2 <- Part$Indices[d2]
    r1 <- A[, 1]; r2 <- A[, 2]

    A[, 4] <- g1; A[, 5] <- g2
    if (id > 1) A <- A[!(g1 %in% Reference_Index & g2 %in% Reference_Index), , drop = FALSE]
    All[[id]] <- A
    if (id == 1) First <- Part$Anchor
  }
  Anchors <- do.call(rbind, All)

  Forward <- Anchors[Anchors[, 4] < Anchors[, 5], , drop = FALSE]
  Forward <- Forward[order(Forward[, 5], Forward[, 4]), , drop = FALSE]
  Mirrored <- Forward; Mirrored[, c(1,2,4,5)] <- Forward[, c(2,1,5,4)]
  Anchors <- rbind(Forward, Mirrored)
  First@anchors <- Anchors
  First@offsets <- Offset
  First@reference.objects <- Reference_Index
  First@object.list <- lapply(Objects, function(x) suppressWarnings(DietSeurat(x, features = Features[1:2])))
  First@command@params$reference <- Reference_Index
  methods::validObject(First)
  say(paste("anchor chunks combined", length(Chunks), "rows", nrow(Anchors)))
  First
}
