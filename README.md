## Pan-fibrosis

This repository contains analysis code for a cross-organ fibroblast atlas of fibrotic diseases. The public scRNA-seq and snRNA-seq atlas includes 4,755,536 cells and nuclei from five organs. Five fibroblast subtypes are described: ACSM3-positive myofibroblasts, inflammatory fibroblasts (infFibs), matrix fibroblasts (matFibs), ACTA2-positive myofibroblasts, and antigen-presenting fibroblasts (apFibs). The four-category framework groups infFibs and matFibs as ECM-associated fibroblasts (ecmFibs).

The analyses characterize transcriptional states, MYC-associated programs, inferred trajectory relationships, predicted ligand–receptor interactions, spatial associations, and clinical associations in the evaluated cohorts. Trajectory estimates do not establish lineage ancestry. Complementary chromatin, spatial, and clinical analyses cover lung, heart, and kidney; they do not establish mechanisms or therapeutic efficacy across all five organs. Spatial compositional clusters were defined separately by organ and are not assumed to be identical cross-organ ecotypes.

### Analysis code

Scripts 01–28 cover the analyses presented in the manuscript and its supplementary material. Files 15 and 16 contain Shell commands and use the `.sh` extension.

| File | Analysis |
| --- | --- |
| 18.Metadata_and_Representation.R | Current subtype nomenclature and dataset, organ, and sequencing-modality composition |
| 19.Leave_One_Organ_Out.R | Four-organ integrations and independent single-organ clustering; marker-based manual annotation |
| 20.Alternative_Integration.R | Executed CCA, RPCA, scVI, LIGER, and Scanorama implementations and supporting code |
| 21–24.Cell_Proportion_*.R | Donor-level, organ-stratified beta-binomial regressions for major cell types, four/five fibroblast categories, and endothelial subtypes |
| 25.Cell_Proportion_Forest.R | Forest plots for the cell-proportion comparisons |
| 26.Clinical_Association_Adjustment.R | Clinical covariate- and fibrosis-adjusted Cox/HC3 linear models with proportional-hazards diagnostics |
| 27.Published_Atlas_Comparison.R | Exploratory overlap comparisons with published fibroblast marker sets |
| 28.Integration_Annotation_Plots.R | Read the existing manual annotations and export method-specific UMAP/marker plots |

### Inputs and execution

R scripts are organized for interactive, block-by-block execution. They require the annotated Seurat objects, sample metadata, gene lists and reference resources named in the scripts. The marker-set workbooks for script 27 must be supplied separately under `Inputs/Atlas_Comparison`; local Excel inputs are excluded from the public repository. Large expression matrices, fragment/loom files and intermediate fitted models are not bundled in this code repository. All paths are relative to the repository. Start R in the repository root before the first analysis; Repo_Dir is retained between interactive blocks. Provide the required input files under data/ as described in DATA_LAYOUT.md. Activate the appropriate Python environment before running Python-backed methods. Generated outputs are stored under Results/ or method-specific checkpoint folders.

Subtype labels were assigned by marker-based annotation. Clustering alone does not reproduce manual annotation decisions. In leave-one-organ-out analyses, held-out organs were independently clustered and annotated; these results are not Anchor-based label-transfer predictions. Alternative methods can retain different cell sets after removal of non-fibroblast clusters.

For cell-proportion inference, each donor contributes the target-cell count and remaining-cell count. Sample_Type = Normal defines the normal group; all other non-empty values define the fibrotic group. Analyses are stratified by organ, restrict datasets to at least two donors per group, include dataset fixed effects where applicable, and use beta-binomial logit models with Wald intervals/tests and Benjamini–Hochberg correction within each of four analysis families. This differs from exploratory pooled cell-level scores or section/spot-level spatial comparisons.

High-complexity QC and marker filtering were used; they are not equivalent to applying a dedicated doublet-detection or ambient-RNA-removal algorithm.

Use the original analysis environments and the session/version records associated with the completed results. The integration support scripts retain the executed implementations and call python from the activated environment.

### How to cite

> “Pan-fibrosis”
