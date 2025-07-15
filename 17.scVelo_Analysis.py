

import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
scv.logging.print_version()

IPF01 = anndata.read_loom(".../Loom/IPF_01.loom")
IPF02 = anndata.read_loom(".../Loom/IPF_02.loom")
IPF03 = anndata.read_loom(".../Loom/IPF_03.loom")
IPF04 = anndata.read_loom(".../Loom/IPF_04.loom")
IPF05 = anndata.read_loom(".../Loom/IPF_05.loom")
HF01 = anndata.read_loom(".../Loom/HF_01.loom")
HF02 = anndata.read_loom(".../Loom/HF_02.loom")
IPF01.var_names_make_unique()
IPF02.var_names_make_unique()
IPF03.var_names_make_unique()
IPF04.var_names_make_unique()
IPF05.var_names_make_unique()
HF01.var_names_make_unique()
HF02.var_names_make_unique()

sample_obs = pd.read_csv(".../scVelo/Input/cellID_obs.csv")
umap = pd.read_csv(".../scVelo/Input/cell_embeddings.csv")
cell_clusters = pd.read_csv(".../scVelo/Input/clusters_obs.csv")

IPF01 = IPF01[np.isin(IPF01.obs.index,sample_obs["x"])]
IPF02 = IPF02[np.isin(IPF02.obs.index,sample_obs["x"])]
IPF03 = IPF03[np.isin(IPF03.obs.index,sample_obs["x"])]
IPF04 = IPF04[np.isin(IPF04.obs.index,sample_obs["x"])]
IPF05 = IPF05[np.isin(IPF05.obs.index,sample_obs["x"])]
HF01 = HF01[np.isin(HF01.obs.index,sample_obs["x"])]
HF02 = HF02[np.isin(HF02.obs.index,sample_obs["x"])]
adata = IPF01.concatenate(IPF02, IPF03, IPF04, IPF05, HF01, HF02)

adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {0:'Cell ID'})
adata_index = adata_index.rename(columns = {"CellID":'Cell ID'})
rep=lambda x : x.split("-")[0]
adata_index["Cell ID"]=adata_index["Cell ID"].apply(rep)

umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'})
umap = umap[np.isin(umap["Cell ID"],adata_index["Cell ID"])]
umap = umap.drop_duplicates(subset = ["Cell ID"])
umap_ordered = adata_index.merge(umap, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
adata.obsm['X_umap'] = umap_ordered.values

cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'Cell ID'})
cell_clusters = cell_clusters[np.isin(cell_clusters["Cell ID"],adata_index["Cell ID"])]
cell_clusters = cell_clusters.drop_duplicates(subset = ["Cell ID"])
cell_clusters_ordered = adata_index.merge(cell_clusters, on = "Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:,1:]
adata.obs['clusters']= np.ravel(cell_clusters_ordered.values)

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap",palette = ["#D6C6A5","#4975A5","#93A95D"],save='.../scVelo/Velocity_Dynamical.svg')
