# Input layout

Run the R scripts from the repository root. Paths are relative; no author-specific drive letters or server directories are required. `Repo_Dir` retains the checkout root while individual scripts change into output folders. If switching checkouts in the same R session, remove `Repo_Dir` and start from the new checkout root.

The following large or external inputs are not included. Use the actual retained data and annotations from the study; do not substitute invented metadata or annotation labels.

| Location | Required input |
| --- | --- |
| data/fibroblast/Fibroblast_Gene_Del_Cluster_Del.RDS | Final annotated fibroblast Seurat object |
| data/organs/ | Lung/Liver/Kidney/Heart/Skin_JoinLayers_Celltype.RDS objects used to prepare leave-one-organ-out inputs |
| data/metadata/sample_info.tsv | Sample_ID, Dataset_ID, Dataset_Type and source sample information; original GB18030 text encoding |
| data/metadata/Meta_Organ.txt | Major-cell metadata, replacing Organ with Lung, Liver, Kidney, Heart or Skin |
| data/metadata/fibroblast_four_subtypes.tsv | Four-category cell metadata |
| data/metadata/fibroblast_five_subtypes.tsv | Five-subtype cell metadata |
| data/metadata/all_subtypes.tsv | Combined cell-type and subtype metadata |
| data/gene_sets/ | Gene_Del.txt and ECM_Gene.txt |
| data/integration/Fibroblast_Preprocessed.RDS | Actual shared input for the alternative-integration comparison |
| data/integration/2.CCA through 7.Scanorama | Completed, manually annotated method-specific objects; used by script 28 |
| data/clinical/ | Fib_Marker50.txt and cohorts/, retaining per-organ and per-dataset subfolders and additional_clinical_information.rds |
| Inputs/Atlas_Comparison/ | Required marker-set workbooks for the published-atlas comparison; not bundled |

Legacy example files under data/ (references, matrices, blacklists and spatial metadata) keep the filenames used in scripts 01–17. Check each script's input calls before running it. Large raw inputs, trained models and generated results are intentionally absent from the code repository. Clinical and spatial metadata must preserve the biological sampling structure described in the manuscript.

Results/ and support checkpoints are generated locally. Run 21–24 before 25 to create the beta-binomial result tables. Script 17 uses the pooled velocity inputs named in that script.
