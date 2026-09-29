import scanpy as sc
import anndata
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import os

print('Loading GSE275199...')
adata = sc.read_10x_mtx(
    'data_raw/external_dtp/gse275199',
    prefix='GSE275199_',
    var_names='gene_symbols',
    cache=False
)
adata.var_names_make_unique()

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find HVGs
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
hvg_genes = adata.var_names[adata.var['highly_variable']]
print(f'HVGs: {len(hvg_genes)}')

# Run PCA
pca = PCA(n_components=20, random_state=42)
X_pca = pca.fit_transform(adata[:, hvg_genes].X.toarray())
adata.obsm['X_pca'] = X_pca
print('Explained variance ratio top 5:', np.round(pca.explained_variance_ratio_[:5], 4))

# Cluster into major groups using KMeans (k=5 to match the 5 known conditions!)
kmeans5 = KMeans(n_clusters=5, random_state=42, n_init=10).fit(X_pca)
adata.obs['kmeans_5'] = [f'Cluster_{c}' for c in kmeans5.labels_]

# Also check k=2 (PC9 vs H358)
kmeans2 = KMeans(n_clusters=2, random_state=42, n_init=10).fit(X_pca[:, :5])
adata.obs['cell_line_kmeans'] = [f'Line_{c}' for c in kmeans2.labels_]

# Markers
test_markers = [
    'EGFR', 'KRAS', 'HOPX', 'NKX2-1', 'SFTPB', 'EPCAM', 'CDH1', 'VIM', 'AXL', 
    'MKI67', 'TOP2A', 'PCNA', 'CDKN1A', 'MT2A', 'YAP1', 'WWTR1', 'DUSP6', 'SPRY4'
]
avail = [m for m in test_markers if m in adata.var_names]

profile5 = pd.DataFrame(index=sorted(adata.obs['kmeans_5'].unique()))
for m in avail:
    idx = adata.var_names.get_loc(m)
    val = np.array(adata.X[:, idx].toarray()).flatten()
    adata.obs[m] = val
    profile5[m] = adata.obs.groupby('kmeans_5', observed=False)[m].mean()

profile5['cell_count'] = adata.obs['kmeans_5'].value_counts()
profile5['median_counts'] = adata.obs.groupby('kmeans_5', observed=False)['total_counts'].median()
profile5['median_genes'] = adata.obs.groupby('kmeans_5', observed=False)['n_genes_by_counts'].median()
profile5['pct_mki67_pos'] = adata.obs.groupby('kmeans_5', observed=False)['MKI67'].apply(lambda x: (x > 0).mean() * 100)

print('\n--- 5-Cluster Profiles ---')
print(np.round(profile5[['cell_count', 'NKX2-1', 'EGFR', 'HOPX', 'MKI67', 'pct_mki67_pos', 'CDKN1A', 'DUSP6', 'SPRY4', 'MT2A']], 3).to_string())

# Line 2 Profiles
profile2 = pd.DataFrame(index=sorted(adata.obs['cell_line_kmeans'].unique()))
for m in avail:
    profile2[m] = adata.obs.groupby('cell_line_kmeans', observed=False)[m].mean()
profile2['cell_count'] = adata.obs['cell_line_kmeans'].value_counts()
print('\n--- 2-Line Profiles (PC9 vs H358) ---')
print(np.round(profile2[['cell_count', 'NKX2-1', 'EGFR', 'HOPX', 'MKI67']], 3).to_string())

os.makedirs('results/external_dtp', exist_ok=True)
profile5.to_csv('results/external_dtp/gse275199_kmeans5_profiles.csv')
adata.obs.to_csv('results/external_dtp/gse275199_obs_with_clusters.csv')
print('\nSaved cluster profiles and obs metadata.')
