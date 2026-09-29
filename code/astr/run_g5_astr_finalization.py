# code/astr/run_g5_astr_finalization.py
# Gate G5: ASTR Finalization, Corrected Terminology, and Multiple-Testing Correction

import os
import json
import anndata
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import mannwhitneyu

os.makedirs('tables/astr', exist_ok=True)
os.makedirs('results/astr', exist_ok=True)

print("--- Running G5: ASTR Finalization across All Conditions ---")
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

# 1. Feature selection: Top 1000 HVGs
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[:, top_gene_indices].toarray()
lib_size = sub_X.sum(axis=1, keepdims=True)
lib_size[lib_size == 0] = 1.0
norm_X = np.log1p((sub_X / lib_size) * 10000.0)

# 2. Strict QC filtering (Corrected Terminology: QC-pass vs QC-excluded)
obs['is_qc_excluded'] = (obs['pct_mito'] >= 15.0) | (obs['total_counts'] < 800) | (obs['n_genes'] < 300)
obs['is_qc_pass'] = ~obs['is_qc_excluded']

qc_pass_obs = obs[obs['is_qc_pass']].copy()

# 3. Compute condition-level dispersion in 5-PC space
condition_records = []
dmso_disp_map = {}

# First pass: compute DMSO baselines
for (cname, plate), grp in qc_pass_obs[qc_pass_obs['drug'] == 'DMSO_TF'].groupby(['cell_name', 'plate'], observed=True):
    n_qc = len(grp)
    if n_qc >= 15:
        pca = PCA(n_components=min(5, n_qc - 1), random_state=42)
        pcs = pca.fit_transform(norm_X[grp['cell_idx']])
        centroid = np.median(pcs, axis=0, keepdims=True)
        dists = np.linalg.norm(pcs - centroid, axis=1)
        dmso_disp_map[(cname, plate)] = {
            'mean_disp': float(np.mean(dists)),
            'dists': dists,
            'th_5pct': float(np.percentile(dists, 95.0))
        }

# Second pass: compute treated conditions
for (cname, drug, plate), grp in qc_pass_obs[qc_pass_obs['drug'] != 'DMSO_TF'].groupby(['cell_name', 'drug', 'plate'], observed=True):
    n_qc = len(grp)
    if n_qc < 25:
        continue
    
    n_ex = len(obs[(obs['cell_name'] == cname) & (obs['plate'] == plate) & (obs['drug'] == drug) & (obs['is_qc_excluded'])])
    n_tot = n_qc + n_ex
    
    pca = PCA(n_components=min(5, n_qc - 1), random_state=42)
    pcs = pca.fit_transform(norm_X[grp['cell_idx']])
    centroid = np.median(pcs, axis=0, keepdims=True)
    dists = np.linalg.norm(pcs - centroid, axis=1)
    
    mean_disp = float(np.mean(dists))
    th_5pct = float(np.percentile(dists, 95.0))
    ast_n = int(np.sum(dists >= th_5pct))
    ast_fraction = float(ast_n / n_qc)
    
    # Matched DMSO baseline
    dmso_info = dmso_disp_map.get((cname, plate), None)
    if dmso_info is not None:
        base_disp = dmso_info['mean_disp']
        dmso_dists = dmso_info['dists']
        delta_disp = mean_disp - base_disp
        disp_ratio = mean_disp / (base_disp + 1e-4)
        
        # Test treated dists vs DMSO dists
        stat, pval = mannwhitneyu(dists, dmso_dists, alternative='greater')
        
        # Matched DMSO tail fraction using treated threshold
        dmso_tail_fraction = float(np.mean(dmso_dists >= th_5pct))
    else:
        base_disp = np.nan
        delta_disp = np.nan
        disp_ratio = np.nan
        pval = np.nan
        dmso_tail_fraction = np.nan
        
    astr_score = float(delta_disp / 2.0) if not np.isnan(delta_disp) else np.nan
    delta_ast = ast_fraction - dmso_tail_fraction if not np.isnan(dmso_tail_fraction) else np.nan
    
    condition_records.append({
        'condition_id': f"{plate}_{cname}_{drug}".replace(' ', '_'),
        'cell_name': cname,
        'drug': drug,
        'plate': plate,
        'n_total': n_tot,
        'n_qc_pass': n_qc,
        'n_qc_excluded': n_ex,
        'pct_qc_excluded': round((n_ex / n_tot) * 100.0, 2),
        'mean_dispersion': round(mean_disp, 4),
        'dmso_baseline_dispersion': round(base_disp, 4) if not np.isnan(base_disp) else np.nan,
        'delta_dispersion': round(delta_disp, 4) if not np.isnan(delta_disp) else np.nan,
        'dispersion_ratio': round(disp_ratio, 4) if not np.isnan(disp_ratio) else np.nan,
        'astr_score': round(astr_score, 4) if not np.isnan(astr_score) else np.nan,
        'ast_cell_count': ast_n,
        'ast_fraction': round(ast_fraction, 4),
        'dmso_tail_fraction': round(dmso_tail_fraction, 4) if not np.isnan(dmso_tail_fraction) else np.nan,
        'delta_ast_fraction': round(delta_ast, 4) if not np.isnan(delta_ast) else np.nan,
        'p_value_dispersion': pval
    })

astr_df = pd.DataFrame(condition_records)

# 4. Benjamini-Hochberg FDR correction
valid_p = astr_df['p_value_dispersion'].dropna()
order = np.argsort(valid_p)
ranks = np.empty_like(order)
ranks[order] = np.arange(1, len(valid_p) + 1)
fdr = valid_p.values * len(valid_p) / ranks
# Monotonicity adjustment
fdr_sorted = np.minimum.accumulate(fdr[order][::-1])[::-1]
fdr_adjusted = np.empty_like(fdr_sorted)
fdr_adjusted[order] = np.clip(fdr_sorted, 0.0, 1.0)

astr_df.loc[valid_p.index, 'fdr_dispersion'] = fdr_adjusted

# 5. Operational significance criteria: delta_disp > 0.50, N_qc >= 30, fdr < 0.05
astr_df['is_significant_tail'] = (
    (astr_df['delta_dispersion'] > 0.50) & 
    (astr_df['n_qc_pass'] >= 30) & 
    (astr_df['fdr_dispersion'] < 0.05)
)

astr_df.to_csv('tables/astr/astr_all_conditions.csv', index=False)
print(f"Generated tables/astr/astr_all_conditions.csv with {len(astr_df)} conditions.")
print(f"Conditions with FDR < 0.05: {(astr_df['fdr_dispersion'] < 0.05).sum()}")
print(f"Conditions with operational threshold (DeltaDisp > 0.50, N >= 30, FDR < 0.05): {astr_df['is_significant_tail'].sum()}")
