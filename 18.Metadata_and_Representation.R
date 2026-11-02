Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

dir.create(file.path(Repo_Dir, "Results/Representation"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(Repo_Dir, "Results/Representation"))
library(Seurat)
library(tidyverse)
Fibroblast <- readRDS(file.path(Repo_Dir, "data/fibroblast/Fibroblast_Gene_Del_Cluster_Del.RDS"))
Sample_Info <- read.delim(file.path(Repo_Dir, "data/metadata/sample_info.tsv"), fileEncoding = "GB18030", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(Sample_Info$Sample_ID))
Sample_Match <- match(Fibroblast$orig.ident, Sample_Info$Sample_ID)
stopifnot(!anyNA(Sample_Match))
Fibroblast$Dataset_ID <- Sample_Info$Dataset_ID[Sample_Match]
Fibroblast$Dataset_Type <- Sample_Info$Dataset_Type[Sample_Match]
if (any(as.character(Fibroblast$Subtype1) == "myoFib")) Fibroblast$Subtype1 <- recode(as.character(Fibroblast$Subtype1), "matFib" = "ACSM3+myoFib", "myoFib" = "ACTA2+myoFib")
if (any(as.character(Fibroblast$Subtype2) %in% c("myoFib", "ecmFib_1", "ecmFib_2"))) Fibroblast$Subtype2 <- recode(as.character(Fibroblast$Subtype2), "matFib" = "ACSM3+myoFib", "myoFib" = "ACTA2+myoFib", "ecmFib_1" = "matFib", "ecmFib_2" = "infFib")
Fibroblast$Subtype1 <- gsub("\\+ +", "+", as.character(Fibroblast$Subtype1))
Fibroblast$Subtype2 <- gsub("\\+ +", "+", as.character(Fibroblast$Subtype2))
Order <- c("ACSM3+myoFib", "infFib", "matFib", "ACTA2+myoFib", "apFib")
stopifnot(all(Fibroblast$Subtype2 %in% Order), all(Fibroblast$Dataset_Type %in% c("scRNA-Seq", "snRNA-Seq")))
Fibroblast$Subtype2 <- factor(Fibroblast$Subtype2, levels = Order)
Tissue_Color <- c("Lung" = "#DCC3D0", "Liver" = "#AB82FF", "Heart" = "#FCBA75", "Kidney" = "#499994", "Skin" = "#E56464")
Dataset_Type_Color <- c("scRNA-Seq" = "#1B9E77", "snRNA-Seq" = "#7570B3")
Dataset_Color <- c(
  "GSE115469" = "#E64B35", "GSE122960" = "#4DBBD5", "GSE128033" = "#00A087",
  "GSE128169" = "#3C5488", "GSE131685" = "#F39B7F", "GSE131882" = "#8491B4",
  "GSE132771" = "#91D1C2", "GSE135893" = "#DC0000", "GSE136103" = "#7E6148",
  "GSE136831" = "#B09C85", "GSE138669" = "#00468B", "GSE140989" = "#ED0000",
  "GSE145154" = "#42B540", "GSE151302" = "#0099B4", "GSE151671" = "#925E9F",
  "GSE156310" = "#FDAF91", "GSE156326" = "#AD002A", "GSE159354" = "#ADB6B6",
  "GSE161470" = "#1B1919", "GSE161685" = "#3B4992", "GSE163973" = "#EE0000",
  "GSE171314" = "#008B45", "GSE174219" = "#631879", "GSE174748" = "#008280",
  "GSE181297" = "#BB0021", "GSE181316" = "#5F559B", "GSE182159" = "#A20056",
  "GSE183276" = "#808180", "GSE183277" = "#0073C2", "GSE185477" = "#EFC000",
  "GSE185778" = "#868686", "GSE186343" = "#CD534C", "GSE189539" = "#7AA6DC",
  "GSE189600" = "#003C67", "GSE195460" = "#8F7700", "GSE195718" = "#3B3B3B",
  "GSE199711" = "#A73030", "GSE209635" = "#4A6990", "GSE211785" = "#BC3C29",
  "GSE212046" = "#0072B5", "GSE212837" = "#E18727", "GSE213017" = "#20854E",
  "GSE214085" = "#7876B1", "GSE214088" = "#6F99AD", "GSE220251" = "#FFDC91",
  "GSE220300" = "#EE4C97", "GSE227136" = "#800000", "GSE249279" = "#767676",
  "GSE253348" = "#FFA319", "SCP1849" = "#8A9045", "Zenodo6962685" = "#155F83"
)
for (Group in c("Tissue", "Dataset_Type", "Dataset_ID")) {
  Plot_Data <- Fibroblast[[]] %>% count(Subtype2, .data[[Group]], name = "Cells") %>% group_by(Subtype2) %>% mutate(Proportion = Cells / sum(Cells)) %>% ungroup()
  P1 <- ggplot(Plot_Data, aes(Subtype2, Proportion, fill = .data[[Group]])) +
    geom_col(width = 0.75) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Cell proportion", fill = Group) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))
  if (Group == "Tissue") P1 <- P1 + scale_fill_manual(values = Tissue_Color)
  if (Group == "Dataset_Type") P1 <- P1 + scale_fill_manual(values = Dataset_Type_Color)
  if (Group == "Dataset_ID") P1 <- P1 + scale_fill_manual(values = Dataset_Color)
  ggsave(paste0("Representation_", Group, ".pdf"), plot = P1, width = 9, height = 5)
  write_tsv(Plot_Data, paste0("Representation_", Group, ".tsv"))
}
