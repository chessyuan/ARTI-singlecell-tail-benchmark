import os
import json
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

os.makedirs('results/external_dtp', exist_ok=True)
os.makedirs('tables/external_dtp', exist_ok=True)

print('--- Gate G9: Independent DTP Reference Construction ---')
adata = anndata.read_h5ad('data_processed/external_dtp/gse275199_processed_adata.h5ad')

adata_qc = adata[adata.obs['qc_pass']].copy()
print(f'QC-pass cells: {adata_qc.n_obs} / {adata.n_obs}')

h358_dtp = adata_qc[adata_qc.obs['condition'] == 'H358_Day8_Sotorasib_500nM_DTP']
h358_ut = adata_qc[adata_qc.obs['condition'] == 'H358_Day0_Untreated']
print(f'H358 DTP cells: {h358_dtp.n_obs}, Untreated cells: {h358_ut.n_obs}')

def get_pseudobulk_mean(sub_adata):
    counts = sub_adata.X
    if hasattr(counts, 'toarray'):
        counts = counts.toarray()
    lib_size = counts.sum(axis=1, keepdims=True)
    lib_size[lib_size == 0] = 1.0
    cpm = (counts / lib_size) * 1e4
    log_cpm = np.log1p(cpm)
    return np.array(np.mean(log_cpm, axis=0)).flatten()

mean_dtp = get_pseudobulk_mean(h358_dtp)
mean_ut = get_pseudobulk_mean(h358_ut)

log2fc_h358 = (mean_dtp - mean_ut) / np.log(2.0)
gene_names = np.array(adata.var_names)

dtp_expr_rate = np.array((h358_dtp.X > 0).mean(axis=0)).flatten()
ut_expr_rate = np.array((h358_ut.X > 0).mean(axis=0)).flatten()
expressed_mask = (dtp_expr_rate >= 0.05) | (ut_expr_rate >= 0.05)

dtp_de_df = pd.DataFrame({
    'gene': gene_names,
    'mean_dtp_logcpm': np.round(mean_dtp, 4),
    'mean_untreated_logcpm': np.round(mean_ut, 4),
    'log2FC': np.round(log2fc_h358, 4)
})

dtp_de_filtered = dtp_de_df[expressed_mask].sort_values(by='log2FC', ascending=False).copy()

top100_up = dtp_de_filtered.head(100)['gene'].tolist()
top100_down = dtp_de_filtered.tail(100)['gene'].tolist()

generic_stress_genes = set([
    'MT2A', 'MT1E', 'MT1X', 'HSPA1A', 'HSPA1B', 'HSP90AA1', 'HSPB1', 'DDIT3', 'ATF4', 
    'ATF3', 'EGR1', 'FOS', 'JUN', 'JUNB', 'IER2', 'IER3', 'BAX', 'PMAIP1', 'BBC3'
])

top100_up_stress_purged = [g for g in top100_up if g not in generic_stress_genes]

# PC9 Trajectory
pc9_ut = adata_qc[adata_qc.obs['condition'] == 'PC9_Day0_Untreated']
pc9_d4 = adata_qc[adata_qc.obs['condition'] == 'PC9_Day4_Osimertinib_150nM_Acute']
pc9_d11 = adata_qc[adata_qc.obs['condition'] == 'PC9_Day11_Osimertinib_150nM_DTP']

mean_pc9_ut = get_pseudobulk_mean(pc9_ut)
mean_pc9_d4 = get_pseudobulk_mean(pc9_d4)
mean_pc9_d11 = get_pseudobulk_mean(pc9_d11)

log2fc_pc9_d4 = (mean_pc9_d4 - mean_pc9_ut) / np.log(2.0)
log2fc_pc9_d11 = (mean_pc9_d11 - mean_pc9_ut) / np.log(2.0)

pc9_expr_rate = np.array((pc9_d11.X > 0).mean(axis=0)).flatten()
shared_mask = expressed_mask & (pc9_expr_rate >= 0.05)

r_h358_pc9_d11, p_h358_pc9_d11 = spearmanr(log2fc_h358[shared_mask], log2fc_pc9_d11[shared_mask])
r_h358_pc9_d4, p_h358_pc9_d4 = spearmanr(log2fc_h358[shared_mask], log2fc_pc9_d4[shared_mask])

print(f'\nCross-System Concordance (H358 Sotorasib DTP vs PC9 Osimertinib DTP):')
print(f'  Shared genes evaluated: {np.sum(shared_mask)}')
print(f'  Spearman rho with PC9 Day 11 DTP: {r_h358_pc9_d11:.4f} (p = {p_h358_pc9_d11:.2e})')
print(f'  Spearman rho with PC9 Day 4 Acute: {r_h358_pc9_d4:.4f} (p = {p_h358_pc9_d4:.2e})')

dtp_signatures = {
    'reference_system': 'GSE275199 NCI-H358 8-day Sotorasib 500 nM',
    'cell_counts': {
        'h358_dtp_8d': int(h358_dtp.n_obs),
        'h358_untreated': int(h358_ut.n_obs),
        'pc9_dtp_11d': int(pc9_d11.n_obs),
        'pc9_acute_4d': int(pc9_d4.n_obs),
        'pc9_untreated': int(pc9_ut.n_obs)
    },
    'top100_dtp_upregulated': top100_up,
    'top100_dtp_downregulated': top100_down,
    'top_dtp_up_stress_purged': top100_up_stress_purged,
    'cross_system_spearman_rho_pc9_dtp': float(r_h358_pc9_d11),
    'cross_system_spearman_rho_pc9_acute': float(r_h358_pc9_d4)
}

with open('results/external_dtp/dtp_reference_signatures.json', 'w') as f:
    json.dump(dtp_signatures, f, indent=2)

dtp_de_filtered.to_csv('tables/external_dtp/dtp_signature_genes.csv', index=False)
print('Gate G9 complete: Independent DTP reference constructed successfully.')
