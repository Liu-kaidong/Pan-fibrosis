Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

dir.create(file.path(Repo_Dir, "Results/Atlas_Comparison"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(Repo_Dir, "Results/Atlas_Comparison"))

rm(list = setdiff(ls(), "Repo_Dir"))
library(readxl)
library(dplyr)
library(ggplot2)
dir.create("Comparison/Four_Subtypes", recursive = TRUE, showWarnings = FALSE)
Our_Markers <- read_excel(file.path(Repo_Dir, "Inputs/Atlas_Comparison/our_four_subtypes.xlsx"), .name_repair = "minimal")
Reference_Markers <- read_excel(file.path(Repo_Dir, "Inputs/Atlas_Comparison/reference_subtypes.xlsx"), .name_repair = "minimal")
Our_List <- lapply(Our_Markers, function(x) unique(trimws(as.character(x[!is.na(x) & nzchar(trimws(as.character(x)))]))))
Reference_List <- lapply(Reference_Markers, function(x) unique(trimws(as.character(x[!is.na(x) & nzchar(trimws(as.character(x)))]))))
Order <- c("ACSM3+ myoFib", "ecmFib", "ACTA2+ myoFib", "apFib")

Our_List <- Our_List[Order]

Background <- sort(unique(unlist(c(Our_List, Reference_List), use.names = FALSE)))
writeLines(Background, "Comparison/Four_Subtypes/Marker_Union_Background_Four_Subtypes.txt", useBytes = TRUE)

Reference_Info <- tibble(
  Reference_Subtype = names(Reference_List),
  Study = rep(c("Buechler 2021", "Tsukui 2024", "Fang 2025", "Amrute 2024", "Kuppe 2021", "Ramachandran 2019", "Steele 2025"), c(4, 4, 4, 13, 9, 2, 6))
)
Gene_Set_QC <- bind_rows(
  tibble(Source = "This study", Subtype = names(Our_List), N_genes = lengths(Our_List)),
  Reference_Info %>% transmute(Source = Study, Subtype = Reference_Subtype, N_genes = lengths(Reference_List))
)
write.csv(Gene_Set_QC, "Comparison/Four_Subtypes/Gene_Set_Sizes_Four_Subtypes.csv", row.names = FALSE, fileEncoding = "UTF-8")

Comparison <- bind_rows(lapply(names(Our_List), function(Our_Subtype) {
  bind_rows(lapply(names(Reference_List), function(Reference_Subtype) {
    N_our <- length(Our_List[[Our_Subtype]])
    N_reference <- length(Reference_List[[Reference_Subtype]])
    Shared <- intersect(Our_List[[Our_Subtype]], Reference_List[[Reference_Subtype]])
    N_overlap <- length(Shared)
    N_background <- length(Background)
    Fisher <- fisher.test(
      matrix(c(N_overlap, N_our - N_overlap, N_reference - N_overlap, N_background - N_our - N_reference + N_overlap), nrow = 2, byrow = TRUE),
      alternative = "greater", conf.int = FALSE)
    tibble(
      Our_Subtype = Our_Subtype, Reference_Subtype = Reference_Subtype,
      N_our = N_our, N_reference = N_reference, N_overlap = N_overlap,
      N_background = N_background, Background_type = "Union_of_input_Top_genes_exploratory",
      Jaccard = N_overlap / (N_our + N_reference - N_overlap),
      Dice = 2 * N_overlap / (N_our + N_reference),
      Overlap_coefficient = N_overlap / min(N_our, N_reference),
      Our_coverage = N_overlap / N_our, Reference_coverage = N_overlap / N_reference,
      Expected_overlap = N_our * N_reference / N_background,
      Fold_enrichment = N_overlap * N_background / (N_our * N_reference),
      Odds_ratio_Fisher = unname(Fisher$estimate),
      P_hyper = phyper(N_overlap - 1, N_reference, N_background - N_reference, N_our, lower.tail = FALSE),
      P_Fisher_greater = Fisher$p.value,
      Shared_genes = paste(Shared, collapse = ";")
    )
  }))
})) %>% left_join(Reference_Info, by = "Reference_Subtype") %>%
  mutate(P_adj_BH = p.adjust(P_hyper, method = "BH"), BH_below_0.05 = P_adj_BH < 0.05,
         Log2_enrichment = ifelse(N_overlap > 0, log2(Fold_enrichment), NA_real_))
write.csv(Comparison, "Comparison/Four_Subtypes/Marker_Overlap_Enrichment_Four_Subtypes.csv", row.names = FALSE, fileEncoding = "UTF-8")
Comparison %>% filter(BH_below_0.05) %>% select(Our_Subtype, Study, Reference_Subtype, N_overlap, Jaccard, Fold_enrichment, P_adj_BH)

Comparison <- read.csv("Comparison/Four_Subtypes/Marker_Overlap_Enrichment_Four_Subtypes.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
Cell_Type_Color <- c("ACSM3+ myoFib" = "#4975A5", "ecmFib" = "#BC8F8F", "ACTA2+ myoFib" = "#93A95D", "apFib" = "#CD6E72")
Study_Color <- c("Buechler 2021" = "#A7AAB6", "Tsukui 2024" = "#DCC3D0", "Fang 2025" = "#DCC3D0", "Amrute 2024" = "#FCBA75",
                 "Kuppe 2021" = "#499994", "Ramachandran 2019" = "#AB82FF", "Steele 2025" = "#E56464")
Reference_Info <- Comparison %>% distinct(Reference_Subtype, Study, N_reference) %>%
  mutate(Study = factor(Study, levels = unique(Study)), X = row_number() + 0.8 * (as.integer(Study) - 1))
Comparison <- Comparison %>% left_join(Reference_Info %>% select(Reference_Subtype, X), by = "Reference_Subtype") %>%
  mutate(Y = match(Our_Subtype, rev(Order)))
Study_Blocks <- Reference_Info %>% group_by(Study) %>%
  summarise(Xmin = min(X) - 0.45, Xmax = max(X) + 0.45, Center = mean(X), .groups = "drop") %>%
  mutate(Label = sub(" ([0-9]{4})$", "\n\\1", as.character(Study)), Alternate = row_number() %% 2 == 1)
Row_Info <- data.frame(Our_Subtype = rev(Order), Y = seq_along(Order))
P1 <- ggplot(Comparison, aes(X, Y)) +
  geom_rect(data = filter(Study_Blocks, Alternate), aes(xmin = Xmin, xmax = Xmax, ymin = 0.5, ymax = length(Order) + 0.5), inherit.aes = FALSE, fill = "#F4F6F9", color = NA) +
  geom_segment(data = Study_Blocks, aes(x = Xmin, xend = Xmax, y = length(Order) + 0.5, yend = length(Order) + 0.5), inherit.aes = FALSE, color = "#D8DEE6", linewidth = 0.3) +
  geom_segment(data = Study_Blocks, aes(x = Xmin, xend = Xmax, y = 0.5, yend = 0.5), inherit.aes = FALSE, color = "#D8DEE6", linewidth = 0.3) +
  geom_segment(data = data.frame(Y = seq(1.5, length(Order) - 0.5, 1)), aes(x = 0.55, xend = max(Reference_Info$X) + 0.45, y = Y, yend = Y), inherit.aes = FALSE, color = "#E6EAF0", linewidth = 0.25) +
  geom_segment(data = Study_Blocks, aes(x = Xmin + 0.1, xend = Xmax - 0.1, y = length(Order) + 0.86, yend = length(Order) + 0.86, color = Study), inherit.aes = FALSE, linewidth = 2.6, lineend = "butt") +
  geom_text(data = Study_Blocks, aes(x = Center, y = length(Order) + 1.4, label = Label), inherit.aes = FALSE, size = 3.1, lineheight = 1.1, family = "Arial", fontface = "bold", color = "#344054") +
  geom_segment(data = Row_Info, aes(x = -0.4, xend = -0.4, y = Y - 0.2, yend = Y + 0.2, color = Our_Subtype), inherit.aes = FALSE, linewidth = 2.7, lineend = "round") +
  geom_point(data = filter(Comparison, N_overlap > 0), aes(size = Jaccard, fill = Log2_enrichment), shape = 21, color = "#697586", stroke = 0.28) +
  geom_text(data = filter(Comparison, BH_below_0.05), aes(label = "*"), size = 3, vjust = 0.65, color = "#202B3C", family = "Arial") +
  scale_x_continuous(breaks = Reference_Info$X, labels = ifelse(Reference_Info$N_reference < 50, paste0(Reference_Info$Reference_Subtype, " [", Reference_Info$N_reference, "]"), Reference_Info$Reference_Subtype), limits = c(-0.9, max(Reference_Info$X) + 0.65), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(Order), labels = rev(Order), limits = c(0.45, length(Order) + 1.95), expand = c(0, 0)) +
  scale_color_manual(values = c(Cell_Type_Color, Study_Color), guide = "none") +
  scale_size_area(max_size = 4.7, breaks = c(0.02, 0.05, 0.10, 0.15), limits = c(0, max(Comparison$Jaccard)), name = "Marker similarity · Jaccard") +
  scale_fill_gradient2(low = "#4E6F93", mid = "#FAF8F4", high = "#B9655F", midpoint = 0, breaks = c(-1, 0, 1, 2, 3, 4), name = "Overlap enrichment · log2(observed / expected)") +
  guides(size = guide_legend(order = 1, title.position = "top", nrow = 1, override.aes = list(fill = "#DCE2EA", color = "#697586")),
         fill = guide_colorbar(order = 2, title.position = "top", barwidth = grid::unit(44, "mm"), barheight = grid::unit(2.7, "mm"), ticks = FALSE)) +
  labs(x = NULL, y = NULL, title = "Marker correspondence across fibroblast atlases (four subtypes)",
       subtitle = paste0("Exploratory hypergeometric enrichment within a pooled Top-gene candidate background (N = ", format(unique(Comparison$N_background), big.mark = ",", trim = TRUE), ")"),
       caption = paste0("* BH-adjusted P < 0.05 across ", nrow(Comparison), " comparisons. Blank: no shared genes. Reference signatures contain 50 genes unless a smaller count is shown in brackets.")) +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 9.1, color = "#344054", margin = margin(t = 5)),
        axis.text.y = element_text(size = 11.5, face = "bold", color = "#273345", margin = margin(r = 8)), axis.ticks = element_blank(),
        plot.title = element_text(size = 15, face = "bold", color = "#243246", margin = margin(b = 5)), plot.subtitle = element_text(size = 10, color = "#64748B", margin = margin(b = 4)),
        plot.caption = element_text(size = 8.5, color = "#64748B", hjust = 0, margin = margin(t = 5)), plot.title.position = "plot", plot.caption.position = "plot",
        legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 9.5, face = "bold", color = "#344054"), legend.text = element_text(size = 9.2, color = "#475569"),
        legend.spacing.x = grid::unit(8, "mm"), legend.box.spacing = grid::unit(3, "mm"), legend.margin = margin(t = 4, b = 4), plot.margin = margin(8, 8, 6, 8))
grDevices::cairo_pdf("Comparison/Four_Subtypes/Marker_Enrichment_Four_Subtypes.pdf", width = 300 / 25.4, height = 150 / 25.4, family = "Arial")
print(P1)
dev.off()

rm(list = setdiff(ls(), "Repo_Dir"))
library(readxl)
library(dplyr)
library(ggplot2)
dir.create("Comparison/Five_Subtypes", recursive = TRUE, showWarnings = FALSE)
Our_Markers <- read_excel(file.path(Repo_Dir, "Inputs/Atlas_Comparison/our_five_subtypes.xlsx"), .name_repair = "minimal")
Reference_Markers <- read_excel(file.path(Repo_Dir, "Inputs/Atlas_Comparison/reference_subtypes.xlsx"), .name_repair = "minimal")
Our_List <- lapply(Our_Markers, function(x) unique(trimws(as.character(x[!is.na(x) & nzchar(trimws(as.character(x)))]))))
Reference_List <- lapply(Reference_Markers, function(x) unique(trimws(as.character(x[!is.na(x) & nzchar(trimws(as.character(x)))]))))
Order <- c("ACSM3+ myoFib", "ACTA2+ myoFib", "apFib", "infFib", "matFib")

Our_List <- Our_List[Order]

Background <- sort(unique(unlist(c(Our_List, Reference_List), use.names = FALSE)))
writeLines(Background, "Comparison/Five_Subtypes/Marker_Union_Background_Five_Subtypes.txt", useBytes = TRUE)

Reference_Info <- tibble(
  Reference_Subtype = names(Reference_List),
  Study = rep(c("Buechler 2021", "Tsukui 2024", "Fang 2025", "Amrute 2024", "Kuppe 2021", "Ramachandran 2019", "Steele 2025"), c(4, 4, 4, 13, 9, 2, 6))
)
Gene_Set_QC <- bind_rows(
  tibble(Source = "This study", Subtype = names(Our_List), N_genes = lengths(Our_List)),
  Reference_Info %>% transmute(Source = Study, Subtype = Reference_Subtype, N_genes = lengths(Reference_List))
)
write.csv(Gene_Set_QC, "Comparison/Five_Subtypes/Gene_Set_Sizes_Five_Subtypes.csv", row.names = FALSE, fileEncoding = "UTF-8")

Comparison <- bind_rows(lapply(names(Our_List), function(Our_Subtype) {
  bind_rows(lapply(names(Reference_List), function(Reference_Subtype) {
    N_our <- length(Our_List[[Our_Subtype]])
    N_reference <- length(Reference_List[[Reference_Subtype]])
    Shared <- intersect(Our_List[[Our_Subtype]], Reference_List[[Reference_Subtype]])
    N_overlap <- length(Shared)
    N_background <- length(Background)
    Fisher <- fisher.test(
      matrix(c(N_overlap, N_our - N_overlap, N_reference - N_overlap, N_background - N_our - N_reference + N_overlap), nrow = 2, byrow = TRUE),
      alternative = "greater", conf.int = FALSE)
    tibble(
      Our_Subtype = Our_Subtype, Reference_Subtype = Reference_Subtype,
      N_our = N_our, N_reference = N_reference, N_overlap = N_overlap,
      N_background = N_background, Background_type = "Union_of_input_Top_genes_exploratory",
      Jaccard = N_overlap / (N_our + N_reference - N_overlap),
      Dice = 2 * N_overlap / (N_our + N_reference),
      Overlap_coefficient = N_overlap / min(N_our, N_reference),
      Our_coverage = N_overlap / N_our, Reference_coverage = N_overlap / N_reference,
      Expected_overlap = N_our * N_reference / N_background,
      Fold_enrichment = N_overlap * N_background / (N_our * N_reference),
      Odds_ratio_Fisher = unname(Fisher$estimate),
      P_hyper = phyper(N_overlap - 1, N_reference, N_background - N_reference, N_our, lower.tail = FALSE),
      P_Fisher_greater = Fisher$p.value,
      Shared_genes = paste(Shared, collapse = ";")
    )
  }))
})) %>% left_join(Reference_Info, by = "Reference_Subtype") %>%
  mutate(P_adj_BH = p.adjust(P_hyper, method = "BH"), BH_below_0.05 = P_adj_BH < 0.05,
         Log2_enrichment = ifelse(N_overlap > 0, log2(Fold_enrichment), NA_real_))
write.csv(Comparison, "Comparison/Five_Subtypes/Marker_Overlap_Enrichment_Five_Subtypes.csv", row.names = FALSE, fileEncoding = "UTF-8")
Comparison %>% filter(BH_below_0.05) %>% select(Our_Subtype, Study, Reference_Subtype, N_overlap, Jaccard, Fold_enrichment, P_adj_BH)

Comparison <- read.csv("Comparison/Five_Subtypes/Marker_Overlap_Enrichment_Five_Subtypes.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
Cell_Type_Color <- c("ACSM3+ myoFib" = "#4975A5", "infFib" = "#D6C6A5", "matFib" = "#62A883", "ACTA2+ myoFib" = "#93A95D", "apFib" = "#CD6E72")
Study_Color <- c("Buechler 2021" = "#A7AAB6", "Tsukui 2024" = "#DCC3D0", "Fang 2025" = "#DCC3D0", "Amrute 2024" = "#FCBA75",
                 "Kuppe 2021" = "#499994", "Ramachandran 2019" = "#AB82FF", "Steele 2025" = "#E56464")
Reference_Info <- Comparison %>% distinct(Reference_Subtype, Study, N_reference) %>%
  mutate(Study = factor(Study, levels = unique(Study)), X = row_number() + 0.8 * (as.integer(Study) - 1))
Comparison <- Comparison %>% left_join(Reference_Info %>% select(Reference_Subtype, X), by = "Reference_Subtype") %>%
  mutate(Y = match(Our_Subtype, rev(Order)))
Study_Blocks <- Reference_Info %>% group_by(Study) %>%
  summarise(Xmin = min(X) - 0.45, Xmax = max(X) + 0.45, Center = mean(X), .groups = "drop") %>%
  mutate(Label = sub(" ([0-9]{4})$", "\n\\1", as.character(Study)), Alternate = row_number() %% 2 == 1)
Row_Info <- data.frame(Our_Subtype = rev(Order), Y = seq_along(Order))
P1 <- ggplot(Comparison, aes(X, Y)) +
  geom_rect(data = filter(Study_Blocks, Alternate), aes(xmin = Xmin, xmax = Xmax, ymin = 0.5, ymax = length(Order) + 0.5), inherit.aes = FALSE, fill = "#F4F6F9", color = NA) +
  geom_segment(data = Study_Blocks, aes(x = Xmin, xend = Xmax, y = length(Order) + 0.5, yend = length(Order) + 0.5), inherit.aes = FALSE, color = "#D8DEE6", linewidth = 0.3) +
  geom_segment(data = Study_Blocks, aes(x = Xmin, xend = Xmax, y = 0.5, yend = 0.5), inherit.aes = FALSE, color = "#D8DEE6", linewidth = 0.3) +
  geom_segment(data = data.frame(Y = seq(1.5, length(Order) - 0.5, 1)), aes(x = 0.55, xend = max(Reference_Info$X) + 0.45, y = Y, yend = Y), inherit.aes = FALSE, color = "#E6EAF0", linewidth = 0.25) +
  geom_segment(data = Study_Blocks, aes(x = Xmin + 0.1, xend = Xmax - 0.1, y = length(Order) + 0.86, yend = length(Order) + 0.86, color = Study), inherit.aes = FALSE, linewidth = 2.6, lineend = "butt") +
  geom_text(data = Study_Blocks, aes(x = Center, y = length(Order) + 1.4, label = Label), inherit.aes = FALSE, size = 3.1, lineheight = 1.1, family = "Arial", fontface = "bold", color = "#344054") +
  geom_segment(data = Row_Info, aes(x = -0.4, xend = -0.4, y = Y - 0.2, yend = Y + 0.2, color = Our_Subtype), inherit.aes = FALSE, linewidth = 2.7, lineend = "round") +
  geom_point(data = filter(Comparison, N_overlap > 0), aes(size = Jaccard, fill = Log2_enrichment), shape = 21, color = "#697586", stroke = 0.28) +
  geom_text(data = filter(Comparison, BH_below_0.05), aes(label = "*"), size = 3, vjust = 0.65, color = "#202B3C", family = "Arial") +
  scale_x_continuous(breaks = Reference_Info$X, labels = ifelse(Reference_Info$N_reference < 50, paste0(Reference_Info$Reference_Subtype, " [", Reference_Info$N_reference, "]"), Reference_Info$Reference_Subtype), limits = c(-0.9, max(Reference_Info$X) + 0.65), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(Order), labels = rev(Order), limits = c(0.45, length(Order) + 1.95), expand = c(0, 0)) +
  scale_color_manual(values = c(Cell_Type_Color, Study_Color), guide = "none") +
  scale_size_area(max_size = 4.7, breaks = c(0.02, 0.05, 0.10, 0.15), limits = c(0, max(Comparison$Jaccard)), name = "Marker similarity · Jaccard") +
  scale_fill_gradient2(low = "#4E6F93", mid = "#FAF8F4", high = "#B9655F", midpoint = 0, breaks = c(-1, 0, 1, 2, 3, 4), name = "Overlap enrichment · log2(observed / expected)") +
  guides(size = guide_legend(order = 1, title.position = "top", nrow = 1, override.aes = list(fill = "#DCE2EA", color = "#697586")),
         fill = guide_colorbar(order = 2, title.position = "top", barwidth = grid::unit(44, "mm"), barheight = grid::unit(2.7, "mm"), ticks = FALSE)) +
  labs(x = NULL, y = NULL, title = "Marker correspondence across fibroblast atlases (five subtypes)",
       subtitle = paste0("Exploratory hypergeometric enrichment within a pooled Top-gene candidate background (N = ", format(unique(Comparison$N_background), big.mark = ",", trim = TRUE), ")"),
       caption = paste0("* BH-adjusted P < 0.05 across ", nrow(Comparison), " comparisons. Blank: no shared genes. Reference signatures contain 50 genes unless a smaller count is shown in brackets.")) +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 9.1, color = "#344054", margin = margin(t = 5)),
        axis.text.y = element_text(size = 11.5, face = "bold", color = "#273345", margin = margin(r = 8)), axis.ticks = element_blank(),
        plot.title = element_text(size = 15, face = "bold", color = "#243246", margin = margin(b = 5)), plot.subtitle = element_text(size = 10, color = "#64748B", margin = margin(b = 4)),
        plot.caption = element_text(size = 8.5, color = "#64748B", hjust = 0, margin = margin(t = 5)), plot.title.position = "plot", plot.caption.position = "plot",
        legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 9.5, face = "bold", color = "#344054"), legend.text = element_text(size = 9.2, color = "#475569"),
        legend.spacing.x = grid::unit(8, "mm"), legend.box.spacing = grid::unit(3, "mm"), legend.margin = margin(t = 4, b = 4), plot.margin = margin(8, 8, 6, 8))
grDevices::cairo_pdf("Comparison/Five_Subtypes/Marker_Enrichment_Five_Subtypes.pdf", width = 300 / 25.4, height = 150 / 25.4, family = "Arial")
print(P1)
dev.off()
