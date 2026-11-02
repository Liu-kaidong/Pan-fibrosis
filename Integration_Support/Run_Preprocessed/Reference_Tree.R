Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)
reference_sample_tree <- function(Anchor) {
  Sizes <- vapply(Anchor@object.list, ncol, integer(1)); N <- length(Sizes)
  Counts <- table(factor(Anchor@anchors$dataset1, levels = seq_len(N)), factor(Anchor@anchors$dataset2, levels = seq_len(N)))
  Counts <- matrix(as.numeric(Counts), nrow = N)
  Similarity <- (Counts + t(Counts)) / outer(Sizes, Sizes, pmin)
  Similarity[upper.tri(Similarity, diag = TRUE)] <- NA_real_
  Reference <- Anchor@reference.objects
  Tree <- Seurat:::BuildSampleTree(similarity.matrix = Similarity[Reference, Reference])
  Seurat:::AdjustSampleTree(Tree, reference.objects = Reference)
}
