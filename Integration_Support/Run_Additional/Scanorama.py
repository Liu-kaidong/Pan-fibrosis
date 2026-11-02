"""Native all-cell Scanorama. Sparse log-normalized HVGs, no sketch/projection."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits
import scanorama

root = Path(__file__).resolve().parent
obs = pd.read_csv(root / 'Checkpoints/Scanorama_Obs.tsv', sep='\t', index_col=0)
obs.index = obs.index.astype(str)
genes = (root / 'Checkpoints/Scanorama_Genes.txt').read_text().splitlines()
with (root / 'Checkpoints/Scanorama_Data.mtx').open('rb') as f:
    x = mmread(f).T.tocsr()
assert x.shape == (507934, 2000) and obs.index.is_unique
batches = sorted(obs['orig.ident'].unique())
indices = [np.flatnonzero(obs['orig.ident'].values == batch) for batch in batches]
datasets = [x[idx].tocsr() for idx in indices]
print('Scanorama all cells:', x.shape, 'batches:', len(batches), flush=True)
params = dict(dimred=30, knn=20, sigma=15, alpha=0.1, approx=True,
              batch_size=5000, hvg=None, seed=123, sketch=False, verbose=2)
(root / 'Logs/Scanorama_Settings.json').write_text(json.dumps(params, indent=2))
with threadpool_limits(limits=4):
    latent, integrated_genes = scanorama.integrate(datasets, [genes for _ in batches], ds_names=batches, **params)
assert len(integrated_genes) == 2000
out = np.empty((len(obs), 30))
for idx, values in zip(indices, latent):
    out[idx] = values
assert np.isfinite(out).all()
pd.DataFrame(out, index=obs.index, columns=[f'SCANORAMA_{i}' for i in range(1,31)]).to_csv(root / 'Checkpoints/Scanorama_Latent.tsv', sep='\t')
print('Scanorama latent checkpoint saved', flush=True)
