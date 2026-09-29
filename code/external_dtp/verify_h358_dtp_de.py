import os
import scanpy as sc
import anndata
import numpy as np
import pandas as pd

obs_df = pd.read_csv('results/external_dtp/gse275199_obs_with_clusters.csv', index_col=0)
adata = sc.read_10x_mtx(
    'data_raw/external_dtp/gse275199',
    prefix='GSE275199_',
    var_names='gene_symbols',
    cache=False
)
adata.var_names_make_unique()
adata.obs = obs_df

# Filter to H358 clusters (Cluster_0: Day 0, Cluster_4: Day 8)
h358_cells = adata[adata.obs['kmeans_5'].isin(['Cluster_0', 'Cluster_4'])].copy()
sc.pp.normalize_total(h358_cells, target_sum=1e4)
sc.pp.log1p(h358_cells)

sc.tl.rank_genes_groups(h358_cells, groupby='kmeans_5', reference='Cluster_0', groups=['Cluster_4'], method='wilcoxon')
de_df = sc.get.rank_genes_groups_df(h358_cells, group='Cluster_4')

os.makedirs('tables/external_dtp', exist_ok=True)
de_df.to_csv('tables/external_dtp/gse275199_h358_dtp_vs_untreated_de.csv', index=False)
print('Successfully saved tables/external_dtp/gse275199_h358_dtp_vs_untreated_de.csv with', len(de_df), 'genes.')
