
conda activate SCENIC
cd .../SCENIC
tfs=.../allTFs_hg38.txt
feather=.../hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
tbl=.../motifs-v9-nr.hgnc-m0.001-o0.0.tbl
pyscenic grn --num_workers 20 --output Result/ecmFib.tsv --method grnboost2 ecmFib.loom  $tfs
pyscenic ctx Result/ecmFib.tsv $feather --annotations_fname $tbl --expression_mtx_fname ecmFib.loom --mode "dask_multiprocessing" --output Result/ecmFib.csv --num_workers 20 --mask_dropouts
pyscenic aucell ecmFib.loom Result/ecmFib.csv --output Result/Result_pySCENIC.loom --num_workers 20

