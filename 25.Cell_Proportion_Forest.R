Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

dir.create(file.path(Repo_Dir, "Results/Proportions"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(Repo_Dir, "Results/Proportions"))

library(tidyverse)

Tissue_Color <- c("Lung" = "#DCC3D0", "Liver" = "#AB82FF", "Heart" = "#FCBA75", "Kidney" = "#499994", "Skin" = "#E56464")
dir.create("主图森林图", showWarnings = FALSE)

for (Figure in c("1e", "2d", "3b", "5h")) {
  Result <- read_tsv(paste0("Fig.", Figure, "/Fig", Figure, "_Split_BetaBinomial_Results.tsv"), na = c("NA", ""), show_col_types = FALSE)
  Result$Cell_Type <- as.character(Result[[switch(Figure, "1e" = "Celltype", "2d" = "Subtype1", "3b" = "Subtype2", "5h" = "Subtype")]])

  if (Figure %in% c("2d", "3b") && any(Result$Cell_Type %in% c("myoFib", "ecmFib_1", "ecmFib_2"))) {
    Result$Cell_Type <- recode(Result$Cell_Type, "matFib" = "ACSM3+myoFib", "myoFib" = "ACTA2+myoFib", "ecmFib_1" = "matFib", "ecmFib_2" = "infFib")
  }

  Order <- switch(Figure,
    "1e" = c("Adipocyte", "B cell", "Cardiomyocyte", "Endothelial cell", "Epithelial cell", "Fibroblast", "Hepatocyte", "Macrophage", "Mast cell", "Melanocyte", "Mural cell", "Neuronal cell", "T cell"),
    "2d" = c("ACSM3+myoFib", "ecmFib", "ACTA2+myoFib", "apFib"),
    "3b" = c("infFib", "matFib"),
    "5h" = c("ArtECs", "CapECs", "LECs", "VenECs")
  )

  Result <- Result %>%
    filter(Cell_Type %in% Order) %>%
    mutate(Cell_Type = factor(Cell_Type, levels = Order),
           Tissue = factor(Tissue, levels = rev(c("Lung", "Liver", "Heart", "Kidney", "Skin"))),
           Estimable = Status == "OK" & is.finite(OR) & is.finite(Lower95) & is.finite(Upper95) &
             OR > 0 & Lower95 > 0 & Upper95 > 0 & Lower95 <= OR & Upper95 >= OR,
           Label = case_when(!Estimable ~ "NE", is.na(FDR) ~ "NA", FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", TRUE ~ "ns"),
           CI_Left = pmax(Lower95, 0),
           CI_Right = pmin(Upper95, 2),
           Label = case_when(Estimable & OR > 2 ~ paste0("OR=", signif(OR, 3), "; ", Label), TRUE ~ Label))

  print(Result %>% select(Tissue, Cell_Type, N_Dataset, Normal_Donor, Fibrosis_Donor, OR, Lower95, Upper95, FDR, Status, Warning))

  P1 <- ggplot(Result, aes(y = Tissue, color = Tissue)) +
    geom_vline(xintercept = 1, linetype = 2, color = "#777777", linewidth = 0.4) +
    geom_segment(data = filter(Result, Estimable, Upper95 >= 0, Lower95 <= 2),
                 aes(x = CI_Left, xend = CI_Right, yend = Tissue), linewidth = 0.6) +
    geom_segment(data = filter(Result, Estimable, Upper95 > 2),
                 aes(x = 1.88, xend = 2, yend = Tissue), linewidth = 0.6,
                 arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed")) +
    geom_point(data = filter(Result, Estimable, OR >= 0, OR <= 2), aes(x = OR), size = 2.4) +
    geom_text(aes(x = 1.95, label = Label), color = "black", size = 2.8, hjust = 1, nudge_y = -0.28) +
    facet_wrap(~ Cell_Type, ncol = if (Figure == "1e") 4 else length(Order)) +
    scale_color_manual(values = Tissue_Color, guide = "none") +
    scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2), labels = c("0", "0.5", "1", "1.5", "2"), expand = expansion(mult = 0)) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian(xlim = c(0, 2)) +
    labs(x = "Fibrosis vs normal: odds ratio (95% CI)", y = NULL,
         caption = "Stars indicate BH-adjusted P values; 95% CIs are not multiplicity-adjusted. NE, not estimable.\nArrows mark CIs extending above 2; ORs above 2 are printed numerically.\n* P_adj < 0.05; ** < 0.01; *** < 0.001. Blank rows have no comparison in the source.") +
    theme_bw(base_size = 10, base_family = "Arial") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "#F4F4F4", color = NA),
          strip.text = element_text(face = "bold"), axis.text = element_text(color = "black"),
          plot.margin = margin(5.5, 10, 5.5, 5.5), panel.spacing = grid::unit(6, "mm"), plot.caption = element_text(size = 8, hjust = 0))

  if (interactive()) print(P1)
  ggsave(paste0("主图森林图/Fig", Figure, "_BetaBinomial_Forest.pdf"), P1, device = cairo_pdf,
         width = if (Figure == "3b") 7 else 12, height = if (Figure == "1e") 10 else 3.8, units = "in")
}
