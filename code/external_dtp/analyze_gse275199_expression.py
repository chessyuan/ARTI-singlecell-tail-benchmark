import scanpy as sc
import anndata
import numpy as np
import pandas as pd
import os

print('Loading GSE275199 matrix...')
adata = sc.read_10x_mtx(
    'data_raw/external_dtp/gse275199',
    prefix='GSE275199_',
    var_names='gene_symbols',
    cache=False
)
adata.var_names_make_unique()

# Calculate basic QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(f'Cells: {adata.n_obs}, Genes: {adata.n_vars}')
print(f'Median total counts: {np.median(adata.obs["total_counts"]):.1f}')
print(f'Median genes detected: {np.median(adata.obs["n_genes_by_counts"]):.1f}')
print(f'Median pct mito: {np.median(adata.obs["pct_counts_mt"]):.2f}%')

# Normalization and HVG
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

print(f'Highly variable genes: {np.sum(adata.var["highly_variable"])}')

# PCA & Neighbors & Leiden
sc.pp.pca(adata, svd_solver='arpack', use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5)

print('Leiden cluster counts:')
print(adata.obs['leiden'].value_counts())

# Check expression of key cell-type and DTP markers
markers = ['EGFR', 'KRAS', 'HOPX', 'AGER', 'SFTPC', 'SFTPA1', 'EPCAM', 'CDH1', 'VIM', 'FN1', 'AXL', 'MET', 'SOX2', 'NGFR', 'YAP1', 'WWTR1', 'MKI67', 'TOP2A', 'PCNA', 'MT2A']
available_markers = [m for m in markers if m in adata.var_names]

cluster_means = pd.DataFrame(index=adata.obs['leiden'].cat.categories)
for m in available_markers:
    # get dense column
    gene_idx = adata.var_names.get_loc(m)
    val = adata.X[:, gene_idx]
    if hasattr(val, 'toarray'):
        val = val.toarray().flatten()
    else:
        val = np.array(val).flatten()
    adata.obs[m] = val
    cluster_means[m] = adata.obs.groupby('leiden', observed=False)[m].mean()

print('\nCluster Mean Expression of Key Markers:')
print(np.round(cluster_means, 3).to_string())

os.makedirs('results/external_dtp', exist_ok=True)
cluster_means.to_csv('results/external_dtp/gse275199_leiden_marker_means.csv')
