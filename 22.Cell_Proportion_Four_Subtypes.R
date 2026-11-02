Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

dir.create(file.path(Repo_Dir, "Results/Proportions/Fig.2d"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(Repo_Dir, "Results/Proportions/Fig.2d"))

rm(list = setdiff(ls(), "Repo_Dir"))
library(tidyverse)
library(glmmTMB)

Fib_Meta <- data.table::fread(file.path(Repo_Dir, "data/metadata/fibroblast_four_subtypes.tsv"), select = c("orig.ident", "Dataset", "Sample_Type", "Fibrosis_State", "Tissue", "Subtype1"), data.table = FALSE)
Result <- Fib_Meta %>%
  transmute(Tissue = as.character(Tissue), Dataset = as.character(Dataset), Sample = as.character(orig.ident),
            Donor_ID = as.character(orig.ident), Disease = as.character(Sample_Type),
            Type = as.character(Fibrosis_State), Subtype1 = as.character(Subtype1)) %>%
  mutate(across(everything(), str_trim), Type = recode(Type, "Firbrosis" = "Fibrosis")) %>%
  count(Tissue, Dataset, Sample, Donor_ID, Disease, Type, Subtype1, name = "Cell_Count")

Group_Mapping <- Result %>% distinct(Tissue, Dataset, Donor_ID, Sample_Type = Disease, Original_Fibrosis_State = Type) %>%
  mutate(Final_Type = if_else(Sample_Type == "Normal", "Normal", "Fibrosis"), Changed = Original_Fibrosis_State != Final_Type,
         Group_Source = "Sample_Type: Normal -> Normal; all other non-empty values -> Fibrosis", Original_Metadata = file.path(Repo_Dir, "data/metadata/fibroblast_four_subtypes.tsv"))
write_tsv(Group_Mapping, "Fig2d_Group_Mapping.tsv")
Result <- Result %>% mutate(Type = if_else(Disease == "Normal", "Normal", "Fibrosis"))
Sample_Info <- Result %>% distinct(Tissue, Dataset, Sample, Donor_ID, Disease, Type)

Result <- Result %>%
  select(Tissue, Dataset, Sample, Donor_ID, Type, Subtype1, Cell_Count) %>%
  group_by(Tissue) %>% complete(nesting(Dataset, Sample, Donor_ID, Type), Subtype1, fill = list(Cell_Count = 0L)) %>%
  group_by(Tissue, Dataset, Sample, Donor_ID, Type) %>%
  mutate(Total_Count = sum(Cell_Count), Freq = Cell_Count / Total_Count) %>% ungroup()
Dataset_Info <- Sample_Info %>%
  count(Tissue, Dataset, Type, name = "N_Donor") %>%
  group_by(Tissue, Dataset) %>% complete(Type = c("Normal", "Fibrosis"), fill = list(N_Donor = 0L)) %>% ungroup() %>%
  pivot_wider(names_from = Type, values_from = N_Donor, values_fill = 0) %>%
  mutate(Eligible = Normal >= 2 & Fibrosis >= 2, Reason = if_else(Eligible, "Included", "Fewer than 2 donors in at least one group"))
Result_All <- Result
Result <- Result %>%
  semi_join(Dataset_Info %>% filter(Eligible), by = c("Tissue", "Dataset")) %>%
  mutate(Type = factor(Type, levels = c("Normal", "Fibrosis")),
         Tissue = factor(Tissue, levels = c("Lung", "Liver", "Kidney", "Heart", "Skin")),
         Study = interaction(Tissue, Dataset, drop = TRUE))

Comparison <- Result %>% distinct(Tissue, Subtype1)
Model_List <- list()
Result_Label <- list()
for (i in seq_len(nrow(Comparison))) {
  Result_Tem <- Result %>% filter(Tissue == Comparison$Tissue[i], Subtype1 == Comparison$Subtype1[i])
  Result_Tem <- droplevels(Result_Tem)
  message("Analysis: ", paste(Comparison$Tissue[i], Comparison$Subtype1[i], sep = "__"))
  Result_Row <- tibble(Tissue = paste(unique(as.character(Result_Tem$Tissue)), collapse = ";"), Subtype1 = Comparison$Subtype1[i],
                       N_Dataset = n_distinct(Result_Tem$Dataset), N_Study = n_distinct(Result_Tem$Study),
                       Normal_Donor = sum(Result_Tem$Type == "Normal"), Fibrosis_Donor = sum(Result_Tem$Type == "Fibrosis"),
                       Normal_Median = median(Result_Tem$Freq[Result_Tem$Type == "Normal"]),
                       Fibrosis_Median = median(Result_Tem$Freq[Result_Tem$Type == "Fibrosis"]),
                       OR = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_, P_Value = NA_real_,
                       Status = "Not fitted", Optimizer = "nlminb", Warning = "")
  if (n_distinct(Result_Tem$Freq) < 2) {
    Result_Row$Status <- "No variation"
    Result_Label[[i]] <- Result_Row
    next
  }
  Model_Formula <- if (n_distinct(Result_Tem$Study) > 1) {
    cbind(Cell_Count, Total_Count - Cell_Count) ~ Type + Study
  } else {
    cbind(Cell_Count, Total_Count - Cell_Count) ~ Type
  }
  Model_Warnings <- character()
  Model <- tryCatch(withCallingHandlers(
    glmmTMB(Model_Formula, family = betabinomial(link = "logit"), data = Result_Tem),
    warning = function(w) {
      Model_Warnings <<- c(Model_Warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (!inherits(Model, "error") && (!isTRUE(Model$sdr$pdHess) || !isTRUE(Model$fit$convergence == 0))) {
    Result_Row$Optimizer <- "BFGS retry"
    Model <- tryCatch(withCallingHandlers(
      glmmTMB(Model_Formula, family = betabinomial(link = "logit"), data = Result_Tem,
              control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))),
      warning = function(w) {
        Model_Warnings <<- c(Model_Warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
  }
  Result_Row$Warning <- paste(unique(Model_Warnings), collapse = " | ")
  if (inherits(Model, "error")) {
    Result_Row$Status <- paste("Fit error:", conditionMessage(Model))
  } else {
    Model_List[[paste(Comparison$Tissue[i], Comparison$Subtype1[i], sep = "__")]] <- Model
    if (!isTRUE(Model$sdr$pdHess) || !isTRUE(Model$fit$convergence == 0)) {
      Result_Row$Status <- "Convergence/Hessian failure"
    } else {
      Coef <- summary(Model)$coefficients$cond["TypeFibrosis", ]
      CI <- exp(Coef["Estimate"] + c(-1, 1) * qnorm(0.975) * Coef["Std. Error"])
      if (all(is.finite(c(Coef, CI, exp(Coef["Estimate"])))) && Coef["Std. Error"] > 0) {
        Result_Row$OR <- unname(exp(Coef["Estimate"]))
        Result_Row$Lower95 <- unname(CI[1])
        Result_Row$Upper95 <- unname(CI[2])
        Result_Row$P_Value <- unname(Coef["Pr(>|z|)"])
        Result_Row$Status <- "OK"
      } else {
        Result_Row$Status <- "Non-finite coefficient or confidence interval"
      }
    }
  }
  Result_Label[[i]] <- Result_Row
}
Result_Label <- bind_rows(Result_Label)
Result_Label$FDR <- NA_real_

Result_Label$FDR[Result_Label$Status == "OK"] <- p.adjust(Result_Label$P_Value[Result_Label$Status == "OK"], method = "BH", n = nrow(Comparison))
Result_Label <- Result_Label %>%
  mutate(P_Label = case_when(is.na(FDR) ~ "NE", FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", TRUE ~ "ns"))
print(Result_Label)
print(Result_Label %>% filter(Status != "OK" | Warning != ""))

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
pdf(NULL)
Plot_List <- list()
Plot_Size <- list()
for (Tissue_Name in c("Lung", "Liver", "Kidney", "Heart", "Skin")) {
  Plot_Data <- Result %>% filter(Tissue == Tissue_Name)
  if (nrow(Plot_Data) == 0) next
  Plot_Label <- Plot_Data %>%
    group_by(Subtype1) %>% summarise(y = max(Freq) + 0.03, .groups = "drop") %>%
    left_join(Result_Label %>% filter(Tissue == Tissue_Name), by = "Subtype1") %>%
    mutate(Label = P_Label)
  P1 <- ggplot(Plot_Data, aes(Type, Freq)) +
    geom_boxplot(fill = NA, color = "black", linewidth = 0.3, outlier.shape = NA, width = 0.65) +
    geom_point(aes(color = Dataset), position = position_jitter(width = 0.2, height = 0, seed = 123), shape = 16, size = 1.1, alpha = 0.6) +
    geom_text(data = Plot_Label, aes(x = 1.5, y = y, label = Label), inherit.aes = FALSE, size = 2.7, vjust = 0) +
    scale_x_discrete(limits = c("Fibrosis", "Normal"), expand = expansion(add = 0.6)) +
    scale_color_manual(values = Dataset_Color, breaks = sort(unique(Plot_Data$Dataset))) +
    guides(color = guide_legend(ncol = max(2, floor(n_distinct(Plot_Data$Subtype1) * 0.75)), byrow = TRUE)) +
    facet_wrap(~ Subtype1, scales = "free_y", nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.3))) +
    labs(x = NULL, y = "Fibroblast subtype proportion", title = Tissue_Name, color = "Dataset") +
    theme_classic(base_size = 9, base_family = "Helvetica") +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "bottom", legend.text = element_text(size = 7), legend.title = element_text(size = 8),
        legend.key.height = grid::unit(3, "mm"), legend.key.width = grid::unit(3, "mm"),
        panel.spacing.x = grid::unit(2.5, "mm"), plot.margin = margin(3, 3, 3, 3))

  Plot_Grob <- ggplotGrob(P1)
  Axis_Columns <- unique(Plot_Grob$layout$l[grepl("^axis-l", Plot_Grob$layout$name)])
  Plot_Grob$widths[Axis_Columns] <- grid::unit(rep(max(grid::convertWidth(Plot_Grob$widths[Axis_Columns], "in", valueOnly = TRUE)), length(Axis_Columns)), "in")
  Plot_Grob$widths[unique(Plot_Grob$layout$l[grepl("^panel", Plot_Grob$layout$name)])] <- grid::unit(0.90, "in")
  Plot_Grob$heights[unique(Plot_Grob$layout$t[grepl("^panel", Plot_Grob$layout$name)])] <- grid::unit(1.45, "in")
  Plot_Width <- grid::convertWidth(sum(Plot_Grob$widths), "in", valueOnly = TRUE)
  Plot_Height <- grid::convertHeight(sum(Plot_Grob$heights), "in", valueOnly = TRUE)
  Plot_List[[Tissue_Name]] <- Plot_Grob
  Plot_Size[[Tissue_Name]] <- c(Plot_Width, Plot_Height)
}
if (interactive()) grid::grid.draw(Plot_List[["Lung"]])
dev.off()
pdf("Fig2d_Split_Proportions.pdf", width = max(map_dbl(Plot_Size, 1)), height = max(map_dbl(Plot_Size, 2)), onefile = TRUE, useDingbats = FALSE)
for (Plot_Grob in Plot_List) {
  grid::grid.newpage()
  grid::grid.draw(Plot_Grob)
}
dev.off()

write_tsv(Sample_Info, "Fig2d_Split_Sample_Info.tsv")
write_tsv(Dataset_Info, "Fig2d_Split_Dataset_Inclusion.tsv")
write_tsv(Result_All, "Fig2d_Split_All_Donor_Proportions.tsv")
write_tsv(Result, "Fig2d_Split_Analysis_Donor_Proportions.tsv")
write_tsv(Result_Label, "Fig2d_Split_BetaBinomial_Results.tsv")
saveRDS(Model_List, "Fig2d_Split_BetaBinomial_Models.rds")
capture.output(sessionInfo(), file = "Fig2d_Split_SessionInfo.txt")
