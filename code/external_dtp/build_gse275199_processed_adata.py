import os
import json
import scanpy as sc
import anndata
import numpy as np
import pandas as pd

os.makedirs('data_processed/external_dtp', exist_ok=True)
os.makedirs('results/external_dtp', exist_ok=True)
os.makedirs('tables/external_dtp', exist_ok=True)

print('Compiling standardized GSE275199 AnnData...')
obs_df = pd.read_csv('results/external_dtp/gse275199_obs_with_clusters.csv', index_col=0)
adata = sc.read_10x_mtx(
    'data_raw/external_dtp/gse275199',
    prefix='GSE275199_',
    var_names='gene_symbols',
    cache=False
)
adata.var_names_make_unique()

cluster_map = {
    'Cluster_0': {
        'cell_line': 'NCI-H358',
        'treatment': 'Untreated',
        'timepoint': 'Day 0',
        'condition': 'H358_Day0_Untreated',
        'is_dtp': False,
        'hashtag_sequence': 'AGTAAGTTCAGCG',
        'driver': 'KRAS_G12C',
        'drug': 'None',
        'dose': '0'
    },
    'Cluster_4': {
        'cell_line': 'NCI-H358',
        'treatment': 'Sotorasib 500 nM',
        'timepoint': 'Day 8',
        'condition': 'H358_Day8_Sotorasib_500nM_DTP',
        'is_dtp': True,
        'hashtag_sequence': 'AAGTATCGTTTCG',
        'driver': 'KRAS_G12C',
        'drug': 'Sotorasib',
        'dose': '500 nM'
    },
    'Cluster_1': {
        'cell_line': 'PC9',
        'treatment': 'Untreated',
        'timepoint': 'Day 0',
        'condition': 'PC9_Day0_Untreated',
        'is_dtp': False,
        'hashtag_sequence': 'GTCAACTCTTTAG',
        'driver': 'EGFR_ex19del',
        'drug': 'None',
        'dose': '0'
    },
    'Cluster_2': {
        'cell_line': 'PC9',
        'treatment': 'Osimertinib 150 nM',
        'timepoint': 'Day 4',
        'condition': 'PC9_Day4_Osimertinib_150nM_Acute',
        'is_dtp': False,
        'hashtag_sequence': 'TGATGGCCTATTG',
        'driver': 'EGFR_ex19del',
        'drug': 'Osimertinib',
        'dose': '150 nM'
    },
    'Cluster_3': {
        'cell_line': 'PC9',
        'treatment': 'Osimertinib 150 nM',
        'timepoint': 'Day 11',
        'condition': 'PC9_Day11_Osimertinib_150nM_DTP',
        'is_dtp': True,
        'hashtag_sequence': 'TTCCGCCTCTCTT',
        'driver': 'EGFR_ex19del',
        'drug': 'Osimertinib',
        'dose': '150 nM'
    }
}

for col in ['cell_line', 'treatment', 'timepoint', 'condition', 'is_dtp', 'hashtag_sequence', 'driver', 'drug', 'dose']:
    obs_df[col] = obs_df['kmeans_5'].map(lambda c: cluster_map[c][col])

# QC pass filter (10x standard)
obs_df['qc_pass'] = (obs_df['pct_counts_mt'] < 15.0) & (obs_df['total_counts'] >= 1000) & (obs_df['n_genes_by_counts'] >= 500)
adata.obs = obs_df

# Save processed anndata
adata.write_h5ad('data_processed/external_dtp/gse275199_processed_adata.h5ad')
print('Saved data_processed/external_dtp/gse275199_processed_adata.h5ad')

# Summary table by condition
summary = adata.obs.groupby(['cell_line', 'timepoint', 'condition', 'drug', 'dose'], observed=True).agg(
    n_cells=('total_counts', 'count'),
    n_cells_qc_pass=('qc_pass', 'sum'),
    median_counts=('total_counts', 'median'),
    median_genes=('n_genes_by_counts', 'median'),
    median_pct_mt=('pct_counts_mt', 'median'),
    mean_mki67=('MKI67', 'mean'),
    mean_hopx=('HOPX', 'mean')
).reset_index()

print('\n--- GSE275199 Condition Summary ---')
print(summary.to_string())
summary.to_csv('tables/external_dtp/gse275199_condition_summary.csv', index=False)

with open('results/external_dtp/gse275199_audit_summary.json', 'w') as f:
    json.dump({
        'accession': 'GSE275199',
        'pmid': '40352726',
        'pmcid': 'PMC12063144',
        'total_cells': int(adata.n_obs),
        'total_genes': int(adata.n_vars),
        'conditions': summary.to_dict(orient='records')
    }, f, indent=2)

print('Audit summary written successfully.')
