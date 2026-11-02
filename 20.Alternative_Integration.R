Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

setwd(file.path(Repo_Dir, "Integration_Support/Run_Preprocessed"))
dir.create("Checkpoints", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)
dir.create("Logs", showWarnings = FALSE)
if (!file.exists("Checkpoints/Common.rds")) source("Prepare.R", encoding = "UTF-8")
for (Method in c("CCA", "RPCA", "scVI")) {
  Sys.setenv(GM_STAGE = Method)
  source("Analysis.R", encoding = "UTF-8")
  gc()
}
setwd(file.path(Repo_Dir, "Integration_Support/Run_Additional"))
dir.create("Checkpoints", showWarnings = FALSE)
dir.create("Tables", showWarnings = FALSE)
dir.create("Logs", showWarnings = FALSE)
file.copy("../Run_Preprocessed/Checkpoints/Input_Info.rds", "Checkpoints/Input_Info.rds", overwrite = FALSE)
for (Method in c("LIGER", "Scanorama")) {
  Sys.setenv(GM_STAGE = Method)
  source("Analysis.R", encoding = "UTF-8")
  gc()
}
Sys.unsetenv("GM_STAGE")
