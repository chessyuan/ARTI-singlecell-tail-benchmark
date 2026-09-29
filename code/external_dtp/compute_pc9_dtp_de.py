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

# PC9 clusters: Cluster_1 (Day 0), Cluster_2 (Day 4), Cluster_3 (Day 11)
pc9_cells = adata[adata.obs['kmeans_5'].isin(['Cluster_1', 'Cluster_2', 'Cluster_3'])].copy()
sc.pp.normalize_total(pc9_cells, target_sum=1e4)
sc.pp.log1p(pc9_cells)

# Day 4 vs Day 0
sc.tl.rank_genes_groups(pc9_cells, groupby='kmeans_5', reference='Cluster_1', groups=['Cluster_2'], method='wilcoxon')
de_d4 = sc.get.rank_genes_groups_df(pc9_cells, group='Cluster_2')
de_d4.to_csv('tables/external_dtp/gse275199_pc9_day4_vs_day0_de.csv', index=False)

# Day 11 vs Day 0
sc.tl.rank_genes_groups(pc9_cells, groupby='kmeans_5', reference='Cluster_1', groups=['Cluster_3'], method='wilcoxon')
de_d11 = sc.get.rank_genes_groups_df(pc9_cells, group='Cluster_3')
de_d11.to_csv('tables/external_dtp/gse275199_pc9_day11_vs_day0_de.csv', index=False)

print('Saved PC9 Day 4 and Day 11 DE tables successfully.')
