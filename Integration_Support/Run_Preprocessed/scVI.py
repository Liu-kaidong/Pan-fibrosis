"""All-cell sparse-count scVI integration. No subtype labels are used."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
import torch
import scvi

root = Path(__file__).resolve().parent
torch.set_num_threads(4)
torch.set_float32_matmul_precision('high')
scvi.settings.seed = 123
obs = pd.read_csv(root / 'Checkpoints/scVI_Obs.tsv', sep='\t', index_col=0)
obs.index = obs.index.astype(str)
obs['orig.ident'] = obs['orig.ident'].astype(str).astype('category')
genes = (root / 'Checkpoints/scVI_Genes.txt').read_text().splitlines()
with (root / 'Checkpoints/scVI_Counts.mtx').open('rb') as stream:
    counts = scipy.io.mmread(stream).T.tocsr().astype(np.float32)
adata = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=genes))
scvi.model.SCVI.setup_anndata(adata, batch_key='orig.ident')
model = scvi.model.SCVI(adata, n_latent=30, n_layers=2, gene_likelihood='nb')
gpu = torch.cuda.is_available()
(root / 'Logs/scVI_Device.json').write_text(json.dumps({'gpu': gpu, 'device': torch.cuda.get_device_name(0) if gpu else 'CPU', 'torch': torch.__version__, 'scvi': scvi.__version__, 'cells': adata.n_obs, 'genes': adata.n_vars, 'max_epochs': 200, 'early_stopping': True, 'n_epochs_kl_warmup': 50, 'float32_matmul_precision': 'high'}, indent=2))
print('Training all cells:', adata.shape, 'GPU:', gpu, flush=True)
model.train(max_epochs=200, accelerator='gpu' if gpu else 'cpu', devices=1, batch_size=1024, early_stopping=True, early_stopping_patience=30, plan_kwargs={'n_epochs_kl_warmup': 50})
model.save(str(root / 'Checkpoints/scVI_Model'), overwrite=True, save_anndata=False)
for name, values in model.history.items():
    if hasattr(values, 'to_csv'):
        values.to_csv(root / ('Tables/scVI_Training_' + name + '.tsv'), sep='\t')
latent = model.get_latent_representation()
pd.DataFrame(latent, index=obs.index, columns=[f'SCVI_{i}' for i in range(1,31)]).to_csv(root / 'Checkpoints/scVI_Latent.tsv', sep='\t')
print('scVI latent checkpoint saved', flush=True)
