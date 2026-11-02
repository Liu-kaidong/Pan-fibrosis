Repo_Dir <- normalizePath(if (exists("Repo_Dir")) Repo_Dir else ".", winslash = "/", mustWork = TRUE)

dir.create(file.path(Repo_Dir, "Results/Clinical"), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(Repo_Dir, "Results/Clinical"))

rm(list = setdiff(ls(), "Repo_Dir"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(survival)
library(sandwich)
library(gridExtra)
set.seed(123)

Input <- file.path(Repo_Dir, "data/clinical")
GEO <- readRDS(file.path(Input, "cohorts", "additional_clinical_information.rds"))
Order <- c("ACSM3+myoFib", "infFib", "matFib", "ACTA2+myoFib", "apFib")
Cell_Type_Color <- c("ACSM3+myoFib" = "#4975A5", "infFib" = "#D6C6A5", "matFib" = "#62A883", "ACTA2+myoFib" = "#93A95D", "apFib" = "#CD6E72")
Fib_Marker50 <- read.delim(file.path(Input, "Fib_Marker50.txt"), check.names = FALSE)
Fib_Marker50 <- list("ACSM3+myoFib" = Fib_Marker50$matFib, "infFib" = Fib_Marker50$ecmFib2, "matFib" = Fib_Marker50$ecmFib1, "ACTA2+myoFib" = Fib_Marker50$myoFib, "apFib" = Fib_Marker50$apFib)
Fib_Marker50$ECM <- read.delim(file.path(Repo_Dir, "data/gene_sets/ECM_Gene.txt"))$ECM_Score
Fib_Marker50 <- lapply(Fib_Marker50, function(x) unique(x[!is.na(x) & nzchar(x)]))

Get_Num <- function(x) suppressWarnings(as.numeric(sub(" .*", "", sub("^[^:]*: *", "", as.character(x)))))
Read_Clinical <- function(Organ, Dataset) read.delim(file.path(Input, "cohorts", Organ, Dataset, paste0(Dataset, "_Sample_Info.txt")), check.names = FALSE)
Read_Raw <- function(Dataset) bind_rows(lapply(list.files(file.path(Input, "cohorts", "Lung", Dataset), "sampleInfo.csv$", recursive = TRUE, full.names = TRUE), function(x) read.csv(x, check.names = FALSE)))
Score_Signatures <- function(Exp, Gene_Sets) {

  Gene_Index <- lapply(Gene_Sets, function(x) which(rownames(Exp) %in% x))

  if (requireNamespace("GSVA", quietly = TRUE)) return(GSVA::gsva(GSVA::ssgseaParam(Exp, Gene_Sets, alpha = 0.25, normalize = TRUE), verbose = FALSE))
  Result <- vapply(seq_len(ncol(Exp)), function(j) {
    Rank <- as.integer(rank(Exp[, j], ties.method = "average"))
    Position <- integer(length(Rank))
    Position[order(Rank, decreasing = TRUE)] <- seq_along(Rank)
    vapply(Gene_Index, function(x) {
      Weight <- Rank[x]^0.25
      sum(Weight * (length(Rank) - Position[x] + 1)) / sum(Weight) - (length(Rank) * (length(Rank) + 1) / 2 - sum(length(Rank) - Position[x] + 1)) / (length(Rank) - length(x))
    }, numeric(1))
  }, numeric(length(Gene_Sets)))
  dimnames(Result) <- list(names(Gene_Sets), colnames(Exp))
  Result / diff(range(Result))
}

Scores <- list()
Coverage <- list()
for (Dataset in c("GSE70866", "GSE28221", "GSE47460", "GSE137570", "GSE217427", "GSE17800", "GSE130036", "GSE160997")) {
  message("Signature scoring: ", Dataset)
  File <- list.files(file.path(Input, "cohorts"), paste0("^", Dataset, "_(Normalized_Exp|matrix|TPM_Exp)\\.txt$"), recursive = TRUE, full.names = TRUE)

  if (Dataset == "GSE217427") {
    Exp <- read.delim(File, check.names = FALSE)
    HGNC <- read.delim(file.path(Input, "cohorts", "Kidney", Dataset, "hgnc_complete_set_20240520.txt"))
    Exp <- Exp[Exp$Description %in% HGNC$symbol[HGNC$locus_group == "protein-coding gene"], ]
    Exp <- aggregate(Exp[, -1], by = list(Gene = Exp$Description), FUN = mean)
    rownames(Exp) <- Exp$Gene
    Exp <- as.matrix(Exp[, -1])
  } else Exp <- as.matrix(read.delim(File, row.names = 1, check.names = FALSE))
  storage.mode(Exp) <- "double"
  Coverage[[Dataset]] <- data.frame(Dataset, Signature = names(Fib_Marker50), Input_Genes = lengths(Fib_Marker50), Matched_Genes = vapply(Fib_Marker50, function(x) sum(x %in% rownames(Exp)), integer(1)))
  Scores[[Dataset]] <- as.data.frame(t(Score_Signatures(Exp, Fib_Marker50))) %>% tibble::rownames_to_column("Expression_ID")
  rm(Exp)
  gc()
}

Raw <- Read_Raw("GSE70866")
GSE70866 <- Read_Clinical("Lung", "GSE70866") %>%
  left_join(Raw %>% transmute(geoID, Age = Get_Num(age), Sex = factor(Get_Num(gender), levels = 0:1, labels = c("Female", "Male")), GAP = Get_Num(gap)), by = "geoID") %>%
  transmute(Expression_ID = geoID, Patient = rawID, Age, Sex, GAP, Time = surTime / 30, Event = as.numeric(surStatus)) %>% left_join(Scores$GSE70866, by = "Expression_ID")

GSE28221 <- Read_Clinical("Lung", "GSE28221") %>%
  left_join(GEO$GSE28221 %>% transmute(geoID, Age = coalesce(Get_Num(`age of blood draw`), Get_Num(age)), Sex = factor(gender), GEO_Event = coalesce(Get_Num(outcome), Get_Num(`survival status`)), GEO_Time = coalesce(Get_Num(`time to outcome`), Get_Num(`years to outcome`), Get_Num(`time to outcomel`), Get_Num(`time to survival`))), by = "geoID") %>%
  mutate(Expression_ID = geoID, Patient = geoID, Time = GEO_Time * 12, Event = GEO_Event) %>% left_join(Scores$GSE28221, by = "Expression_ID")
Endpoint_Check <- GSE28221 %>% transmute(Patient, Local_Event = surStatus, GEO_Event, Local_Time_Years = surTime, GEO_Time_Years = GEO_Time, Event_Differs = Local_Event != GEO_Event)

Raw <- Read_Raw("GSE47460")
GSE47460 <- Read_Clinical("Lung", "GSE47460") %>%
  left_join(Raw %>% transmute(geoID, Age = Get_Num(age), Sex = factor(sub("Sex: [12]-", "", gender)), ILD_Subtype = sub("ild subtype: *", "", ild_subtype)), by = "geoID") %>%
  filter(ILD_Subtype == "2-UIP/IPF") %>%
  transmute(Expression_ID = geoID, Patient = rawID, Age, Sex, FEV1 = predicted_fev1_pre_bd, FVC = predicted_fvc_pre_bd, DLCO = predicted_dlco) %>% left_join(Scores$GSE47460, by = "Expression_ID")

GSE137570 <- Read_Clinical("Kidney", "GSE137570") %>%
  left_join(GEO$GSE137570 %>% transmute(Sample_geo_accession = geoID, Age = Get_Num(age), Sex = factor(sex), Hist_Fibrosis = Get_Num(`% tubulointerstitial fibrosis on kidney biopsy`)), by = "Sample_geo_accession") %>%
  transmute(Expression_ID = `SAMPLE ID`, Patient = `SAMPLE ID`, Age, Sex, Hist_Fibrosis, eGFR = GFR) %>% left_join(Scores$GSE137570, by = "Expression_ID")

GSE217427 <- Read_Clinical("Kidney", "GSE217427") %>%
  left_join(GEO$GSE217427 %>% transmute(Sample_geo_accession = geoID, Patient = donor, Age = Get_Num(age), Sex = factor(sex), Region = sub("Kidney ", "", tissue)), by = "Sample_geo_accession") %>%
  mutate(Expression_ID = Sample_title, eGFR = egfr) %>% left_join(Scores$GSE217427, by = "Expression_ID") %>%
  group_by(Patient) %>% summarise(across(all_of(c("Age", "eGFR", Order, "ECM")), mean), Sex = first(Sex), Tissue_Samples = n(), .groups = "drop")

Lines <- readLines(file.path(Input, "cohorts", "Heart", "GSE17800", "GSE17800_series_matrix.txt"))
Get_Series_Row <- function(Key) scan(text = sub("^[^\t]*\t", "", Lines[grepl(Key, Lines)]), what = "", sep = "\t", quote = "\"", quiet = TRUE)
Raw <- data.frame(Sample_geo_accession = Get_Series_Row("^!Sample_geo_accession"), Age = Get_Num(Get_Series_Row("^!Sample_characteristics_ch1.*age:")), Sex = factor(sub("gender: ", "", Get_Series_Row("^!Sample_characteristics_ch1.*gender:"))))
GSE17800 <- Read_Clinical("Heart", "GSE17800") %>% left_join(Raw, by = "Sample_geo_accession") %>%
  transmute(Expression_ID = Sample_geo_accession, Patient = Sample_geo_accession, Age, Sex, LVEF = lvef) %>% left_join(Scores$GSE17800, by = "Expression_ID")

Raw <- as.data.frame(read_excel(file.path(Input, "cohorts", "Heart", "GSE130036", "Table S1.xlsx"), skip = 1))
GSE130036 <- Read_Clinical("Heart", "GSE130036") %>% filter(grepl("^HCM", Sample_title)) %>%
  left_join(Raw %>% transmute(Sample_title = Sample, Age = suppressWarnings(as.numeric(Age)), Sex = factor(Sex)), by = "Sample_title") %>%
  transmute(Expression_ID = Sample_geo_accession, Patient = Sample_title, Age, Sex, LVEF) %>% left_join(Scores$GSE130036, by = "Expression_ID")

Raw <- as.data.frame(read_excel(file.path(Input, "cohorts", "Heart", "GSE160997", "Supplementary Data 1.xlsx"), skip = 2, .name_repair = "unique"))
Raw <- Raw %>% transmute(Patient = paste0("HCM", `HCM Patient`), Age = as.numeric(`Age...2`), Sex = factor(`Gender (M=1)`, levels = 0:1, labels = c("Female", "Male")), NYHA = suppressWarnings(as.numeric(`Pre-myectomy NYHA FC`)), MRI_Fibrosis = suppressWarnings(as.numeric(`LGE %`)))
GSE160997 <- Read_Clinical("Heart", "GSE160997") %>% mutate(Patient = sub(".* - ", "", Sample_title)) %>% left_join(Raw, by = "Patient") %>%
  transmute(Expression_ID = Sample_geo_accession, Patient, Age, Sex, NYHA, MRI_Fibrosis, LVEF) %>% left_join(Scores$GSE160997, by = "Expression_ID")

Clinical <- list(GSE70866 = GSE70866, GSE28221 = GSE28221, GSE47460 = GSE47460, GSE137570 = GSE137570,
                 GSE217427 = GSE217427,
                 GSE17800 = GSE17800, GSE130036 = GSE130036, GSE160997 = GSE160997)
Clinical_Result <- list()
Clinical_Input <- list()
Analysis_Check <- list()
Model_Warnings <- character()
Run_Cox <- function(Formula, Data) {
  Unstable <- FALSE
  Fit <- tryCatch(withCallingHandlers({
    coxph(Formula, data = Data, x = TRUE)
  }, warning = function(w) {
    Model_Warnings <<- c(Model_Warnings, conditionMessage(w))
    if (grepl("infinite|converge|iterations", conditionMessage(w))) Unstable <<- TRUE
    invokeRestart("muffleWarning")
  }), error = function(e) NULL)
  if (Unstable || is.null(Fit) || any(!is.finite(coef(Fit)))) NULL else Fit
}
Get_Effect <- function(Fit, Kind) {
  if (Kind == "Survival") {
    S <- summary(Fit)
    return(data.frame(Estimate = exp(coef(Fit)["Score"]), Lower = S$conf.int["Score", "lower .95"], Upper = S$conf.int["Score", "upper .95"], P = S$coefficients["Score", "Pr(>|z|)"]))
  }
  V <- vcovHC(Fit, type = "HC3")
  SE <- sqrt(V["Score", "Score"])
  data.frame(Estimate = coef(Fit)["Score"], Lower = coef(Fit)["Score"] - qt(0.975, Fit$df.residual) * SE, Upper = coef(Fit)["Score"] + qt(0.975, Fit$df.residual) * SE, P = 2 * pt(abs(coef(Fit)["Score"] / SE), Fit$df.residual, lower.tail = FALSE))
}

for (Dataset in names(Clinical)) {
  Endpoints <- if (Dataset %in% c("GSE70866", "GSE28221")) "Survival" else if (Dataset == "GSE47460") c("FEV1", "FVC", "DLCO") else if (grepl("GSE137570|GSE217427", Dataset)) "eGFR" else "LVEF"
  Kind <- if (Endpoints[1] == "Survival") "Survival" else "Function"
  Clinical_Vars <- if (Dataset == "GSE70866") "GAP" else if (Dataset == "GSE160997") c("Age", "Sex", "NYHA") else c("Age", "Sex")
  Fibrosis_Var <- if (Dataset == "GSE137570") "Hist_Fibrosis" else "ECM"
  Models <- list(Unadjusted = character(), AgeSex = c("Age", "Sex"))
  if (Dataset == "GSE70866") Models$GAP <- "GAP"
  if (Dataset == "GSE160997") Models$AgeSex_NYHA <- Clinical_Vars
  Models$FullyAdjusted <- c(Clinical_Vars, Fibrosis_Var)
  for (Endpoint in Endpoints) {
    message("Clinical analysis: ", Dataset, " / ", Endpoint)
    Data <- Clinical[[Dataset]]
    if (Kind == "Function") Data$Outcome <- Data[[Endpoint]]
    Variables <- unique(c("Age", "Sex", Clinical_Vars, Fibrosis_Var, Order, if (Kind == "Survival") c("Time", "Event") else "Outcome"))
    Before <- nrow(Data)
    Data <- Data[complete.cases(Data[, Variables]), ]
    if (Kind == "Survival") Data <- Data[Data$Time > 0 & Data$Event %in% 0:1, ]

    for (Signature in c(Order, "ECM")) Data[[paste0(Signature, "_z")]] <- as.numeric(scale(Data[[Signature]]))
    Clinical_Input[[paste(Dataset, Endpoint)]] <- mutate(Data, Dataset = Dataset, Endpoint = Endpoint)
    Analysis_Check[[paste(Dataset, Endpoint)]] <- data.frame(Dataset, Endpoint, Available_Patients = Before, Analyzed_Patients = nrow(Data), Excluded_Missing = Before - nrow(Data), Events = if (Kind == "Survival") sum(Data$Event) else NA_real_, Clinical_Adjustment = paste(Clinical_Vars, collapse = " + "), Fibrosis_Adjustment = Fibrosis_Var)
    for (Model in names(Models)) {
      Covariates <- Models[[Model]]
      if ("ECM" %in% Covariates) Covariates[Covariates == "ECM"] <- "ECM_z"
      for (Subtype in Order) {
        Data$Score <- Data[[paste0(Subtype, "_z")]]
        Formula <- as.formula(paste(if (Kind == "Survival") "Surv(Time, Event)" else "Outcome", "~", paste(c("Score", Covariates), collapse = " + ")))
        Fit <- if (Kind == "Survival") Run_Cox(Formula, Data) else lm(Formula, data = Data)
        Effect <- Get_Effect(Fit, Kind)
        PH <- if (Kind == "Survival") cox.zph(Fit)$table else NULL
        Clinical_Result[[paste(Dataset, Endpoint, Model, Subtype)]] <- cbind(data.frame(Dataset, Endpoint, Kind, Model, Subtype, N = nrow(Data), Events = if (Kind == "Survival") sum(Data$Event) else NA_real_, Formula = paste(deparse(Formula), collapse = " "), PH_Score_P = if (Kind == "Survival") PH["Score", "p"] else NA_real_, PH_Global_P = if (Kind == "Survival") PH["GLOBAL", "p"] else NA_real_), Effect)
      }
    }

  }
}

Clinical_Result <- bind_rows(Clinical_Result) %>% group_by(Dataset, Endpoint, Model) %>% mutate(P_adj_BH = p.adjust(P, "BH")) %>% ungroup()
Analysis_Check <- bind_rows(Analysis_Check)
Clinical_Input <- bind_rows(Clinical_Input)

Plots <- list()
Adjustment_Panels <- list()
Plot_Data <- Clinical_Result %>% filter(Subtype == "infFib") %>%
  mutate(Model = factor(Model, levels = c("Unadjusted", "AgeSex", "GAP", "AgeSex_NYHA", "FullyAdjusted")))

Plot_InfFib_Adjustment <- function(Data, Title, Subtitle, X_Label, Limits) {
  Data <- Data %>% arrange(Model) %>% mutate(Model_Label = as.character(Model))
  Data$Model_Label[Data$Model == "Unadjusted"] <- "Unadjusted"
  Data$Model_Label[Data$Model == "AgeSex"] <- "Age + sex"
  Data$Model_Label[Data$Model == "GAP"] <- "GAP"
  Data$Model_Label[Data$Model == "AgeSex_NYHA"] <- "Age + sex + NYHA"
  Data$Model_Label[Data$Model == "FullyAdjusted"] <- if (Data$Dataset[1] == "GSE70866") "GAP + ECM" else if (Data$Dataset[1] == "GSE137570") "Age + sex +\nbiopsy fibrosis" else if (Data$Dataset[1] == "GSE160997") "Age + sex +\nNYHA + ECM" else "Age + sex + ECM"
  Data$Model_Label <- factor(Data$Model_Label, levels = rev(Data$Model_Label))
  Data$Display <- ifelse(Data$Model == "FullyAdjusted", "Full", "Other")
  P <- ggplot(Data, aes(Estimate, Model_Label, color = Display)) +
    geom_vline(xintercept = ifelse(Data$Kind[1] == "Survival", 1, 0), linetype = 2, color = "grey65", linewidth = 0.4) +
    geom_errorbar(aes(xmin = Lower, xmax = Upper), orientation = "y", width = 0.12, linewidth = 0.55) +
    geom_point(aes(fill = Display), shape = 21, size = 2.5, stroke = 0.35) +
    scale_color_manual(values = c(Full = "#4975A5", Other = "#8797A5"), guide = "none") +
    scale_fill_manual(values = c(Full = "#4975A5", Other = "#B8C4CD"), guide = "none") +
    scale_y_discrete(expand = expansion(add = 0.55)) +
    labs(title = Title, subtitle = Subtitle, x = X_Label, y = NULL) +
    theme_classic(base_size = 9, base_family = "sans") +
    theme(axis.text = element_text(size = 8, color = "black"), axis.text.y = element_text(lineheight = 0.95, margin = margin(r = 5)),
          axis.title.x = element_text(size = 8, margin = margin(t = 6)), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.x = element_line(linewidth = 0.35), axis.ticks.x = element_line(linewidth = 0.35),
          plot.title = element_text(size = 10, face = "bold", margin = margin(b = 3)), plot.subtitle = element_text(size = 8, margin = margin(b = 8)),
          plot.title.position = "plot", plot.margin = margin(8, 12, 8, 8))
  if (Data$Kind[1] == "Survival") {

    P <- P + scale_x_log10(limits = Limits, breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4))
  } else P <- P + scale_x_continuous(limits = Limits, breaks = scales::breaks_pretty(n = 3))
  P
}

Limits <- range(c(Plot_Data$Lower[Plot_Data$Kind == "Survival"], Plot_Data$Upper[Plot_Data$Kind == "Survival"], 1))
Limits <- exp(log(Limits) + c(-0.08, 0.08) * diff(log(Limits)))
for (Dataset in c("GSE70866", "GSE28221")) {
  Data <- filter(Plot_Data, .data$Dataset == .env$Dataset)
  Adjustment_Panels[[Dataset]] <- Plot_InfFib_Adjustment(Data, paste0(Dataset, ifelse(Dataset == "GSE70866", " | Overall survival", " | Transplant-free survival")),
    paste0("IPF | ", ifelse(Dataset == "GSE70866", "BAL", "PBMC"), " | n = ", Data$N[1], "; events = ", Data$Events[1]), "HR per SD of infFib score (95% CI)", Limits)
}

for (Organ in c("Lung", "Kidney", "Heart")) {
  Datasets <- if (Organ == "Lung") "GSE47460" else if (Organ == "Kidney") c("GSE137570", "GSE217427") else c("GSE17800", "GSE130036", "GSE160997")
  Data <- filter(Plot_Data, Dataset %in% Datasets)
  Limits <- range(c(Data$Lower, Data$Upper, 0))
  Limits <- Limits + c(-0.06, 0.06) * diff(Limits)
  for (Dataset in Datasets) {
    Endpoints <- if (Organ == "Lung") c("FVC", "FEV1", "DLCO") else if (Organ == "Kidney") "eGFR" else "LVEF"
    for (Endpoint in Endpoints) {
      Data <- filter(Plot_Data, .data$Dataset == .env$Dataset, .data$Endpoint == .env$Endpoint)
      Subtitle <- paste0(Dataset, " | n = ", Data$N[1])
      if (Dataset == "GSE47460") Subtitle <- paste0(Subtitle, " | UIP/IPF")
      if (Dataset %in% c("GSE17800", "GSE130036", "GSE160997")) Subtitle <- paste0(Subtitle, ifelse(Dataset == "GSE17800", " | DCM", " | HCM"))
      X_Label <- if (Organ == "Lung") "Difference in % predicted per SD" else if (Organ == "Heart") "LVEF difference (percentage points) per SD" else "eGFR difference per SD"
      Adjustment_Panels[[paste(Dataset, Endpoint, sep = "_")]] <- Plot_InfFib_Adjustment(Data, paste(Organ, Endpoint, sep = " | "), Subtitle, X_Label, Limits)
    }
  }
}

Arrange_Clinical_Panels <- function(Panels, Ncol, Caption, Layout = NULL) {
  Grobs <- lapply(Panels, ggplotGrob)
  Widths <- Reduce(grid::unit.pmax, lapply(Grobs, function(x) x$widths))
  for (i in seq_along(Grobs)) {
    Grobs[[i]]$widths <- Widths
    j <- which(Grobs[[i]]$layout$name == "panel")
    Grobs[[i]]$grobs[[j]] <- grid::grobTree(Grobs[[i]]$grobs[[j]], vp = grid::viewport(name = paste0("clinical-axes-", i)))
  }
  if (is.null(Layout)) Layout <- matrix(seq_along(Grobs), ncol = Ncol, byrow = TRUE)
  gridExtra::arrangeGrob(grobs = Grobs, ncol = Ncol, layout_matrix = Layout, bottom = grid::textGrob(Caption, x = 0, hjust = 0, gp = grid::gpar(fontsize = 8, fontfamily = "sans", lineheight = 1.3)))
}
Plots$Adjustment_Survival <- Arrange_Clinical_Panels(Adjustment_Panels[1:2], 2,
  "Points: HR; bars: 95% CI. Dark blue: fullest available adjustment. Same patients across models.\nGSE28221: constant-HR summaries only; proportional hazards assumption not met. Public GEO TFS endpoint used.")
Plots$Adjustment_Function <- Arrange_Clinical_Panels(Adjustment_Panels[-(1:2)], 3,
  "Points: coefficients per SD of infFib score; bars: 95% HC3 CI. Dark blue: fullest available adjustment.\nSame patients across models. GSE217427: all available tissue samples averaged within each donor (n = 16).", Layout = matrix(c(1:5, NA, 6:8), ncol = 3, byrow = TRUE))

grDevices::pdf("InfFib_Adjustment.pdf", width = 11.7, height = 8.3, family = "Helvetica", useDingbats = FALSE)
grid::grid.newpage()
grid::grid.draw(Plots$Adjustment_Survival)
grid::grid.newpage()
grid::grid.draw(Plots$Adjustment_Function)
dev.off()

dir.create("InfFib_Adjustment_Panels", showWarnings = FALSE)
for (Name in names(Adjustment_Panels)) {
  P <- Adjustment_Panels[[Name]]
  if (Name == "GSE28221") P <- P + labs(caption = "Constant-HR summary; PH assumption not met.") + theme(plot.caption = element_text(size = 7, hjust = 0))
  grDevices::pdf(file.path("InfFib_Adjustment_Panels", paste0(Name, ".pdf")), width = 4.2, height = 2.8, family = "Helvetica", useDingbats = FALSE)
  print(P)
  dev.off()
}

Clinical_Result %>% filter(Subtype == "infFib", Model == "FullyAdjusted") %>% select(Dataset, Endpoint, N, Estimate, Lower, Upper, P_adj_BH)
