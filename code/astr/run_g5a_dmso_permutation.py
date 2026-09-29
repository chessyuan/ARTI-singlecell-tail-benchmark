# code/astr/run_g5a_dmso_permutation.py
# Gate G5A: DMSO Pseudo-Perturbation Permutation Null Calibration

import os
import json
import anndata
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

os.makedirs('results/astr', exist_ok=True)

print("--- Running Gate G5A: DMSO Pseudo-Perturbation Null Calibration ---")
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

# QC filtering
obs['is_qc_excluded'] = (obs['pct_mito'] >= 15.0) | (obs['total_counts'] < 800) | (obs['n_genes'] < 300)
obs['is_qc_pass'] = ~obs['is_dying'] if 'is_dying' in obs.columns else ~obs['is_qc_excluded']

dmso_obs = obs[obs['is_qc_pass'] & (obs['drug'] == 'DMSO_TF')].copy()

# Top 1000 HVGs
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[:, top_gene_indices].toarray()
lib_size = sub_X.sum(axis=1, keepdims=True)
lib_size[lib_size == 0] = 1.0
norm_X = np.log1p((sub_X / lib_size) * 10000.0)

print(f"Total QC-pass DMSO cells available: {len(dmso_obs)}")

np.random.seed(42)
n_permutations = 500
null_records = []

# Cell line groups with sufficient DMSO cells (N >= 60)
eligible_dmso_groups = []
for (cname, plate), grp in dmso_obs.groupby(['cell_name', 'plate'], observed=True):
    if len(grp) >= 60:
        eligible_dmso_groups.append(((cname, plate), grp))

print(f"Eligible DMSO (cell_line, plate) groups for permutation: {len(eligible_dmso_groups)}")

# Run pseudo-perturbations
for perm_i in range(n_permutations):
    # Randomly select an eligible DMSO group
    (cname, plate), grp = eligible_dmso_groups[perm_i % len(eligible_dmso_groups)]
    indices = grp['cell_idx'].values
    n_sample = min(40, len(indices) // 2)
    
    perm_sub = np.random.choice(indices, size=n_sample * 2, replace=False)
    pseudo_treated_idx = perm_sub[:n_sample]
    pseudo_control_idx = perm_sub[n_sample:]
    
    # 5-PC PCA on pseudo-treated
    pca_t = PCA(n_components=min(5, n_sample - 1), random_state=42)
    pcs_t = pca_t.fit_transform(norm_X[pseudo_treated_idx])
    centroid_t = np.median(pcs_t, axis=0, keepdims=True)
    dists_t = np.linalg.norm(pcs_t - centroid_t, axis=1)
    disp_t = float(np.mean(dists_t))
    
    # 5-PC PCA on pseudo-control
    pca_c = PCA(n_components=min(5, n_sample - 1), random_state=42)
    pcs_c = pca_c.fit_transform(norm_X[pseudo_control_idx])
    centroid_c = np.median(pcs_c, axis=0, keepdims=True)
    dists_c = np.linalg.norm(pcs_c - centroid_c, axis=1)
    disp_c = float(np.mean(dists_c))
    
    null_delta_disp = disp_t - disp_c
    null_astr = null_delta_disp / 2.0
    
    # Tail frequency at nominal 95% quantile
    th_5 = np.percentile(dists_t, 95.0)
    null_ast_freq = float(np.mean(dists_t >= th_5))
    
    null_records.append({
        'permutation_id': perm_i + 1,
        'cell_name': cname,
        'plate': plate,
        'n_cells_per_arm': n_sample,
        'pseudo_treated_dispersion': round(disp_t, 4),
        'pseudo_control_dispersion': round(disp_c, 4),
        'null_delta_dispersion': round(null_delta_disp, 4),
        'null_astr_score': round(null_astr, 4),
        'null_ast_fraction': round(null_ast_freq, 4)
    })

null_df = pd.DataFrame(null_records)
null_df.to_csv('results/astr/dmso_null_calibration.csv', index=False)

null_delta = null_df['null_delta_dispersion']
p95_null_delta = float(np.percentile(null_delta, 95.0))
p99_null_delta = float(np.percentile(null_delta, 99.0))
max_null_delta = float(np.max(null_delta))
mean_null_delta = float(np.mean(null_delta))
std_null_delta = float(np.std(null_delta))

# Compare with real conditions from tables/astr/astr_all_conditions.csv
real_df = pd.read_csv('tables/astr/astr_all_conditions.csv')
real_deltas = real_df['delta_dispersion'].dropna()

# Compute empirical p-value for real conditions against DMSO null
real_empirical_p = []
for d in real_df['delta_dispersion']:
    if np.isnan(d):
        real_empirical_p.append(np.nan)
    else:
        emp_p = float(np.mean(null_delta >= d))
        real_empirical_p.append(round(emp_p, 5))

real_df['empirical_p_dmso_null'] = real_empirical_p
real_df.to_csv('tables/astr/astr_all_conditions.csv', index=False)

# False discovery rate of operational threshold (delta_disp > 0.50) in null
null_false_pos_rate_05 = float(np.mean(null_delta > 0.50))

summary_g5a = {
    'n_permutations': n_permutations,
    'mean_null_delta_dispersion': round(mean_null_delta, 4),
    'std_null_delta_dispersion': round(std_null_delta, 4),
    'p95_null_delta_dispersion': round(p95_null_delta, 4),
    'p99_null_delta_dispersion': round(p99_null_delta, 4),
    'max_null_delta_dispersion': round(max_null_delta, 4),
    'operational_threshold_delta_disp': 0.50,
    'null_false_positive_rate_at_050': round(null_false_pos_rate_05, 4),
    'real_conditions_tested': len(real_df),
    'real_conditions_exceeding_p95_null': int(np.sum(real_deltas > p95_null_delta)),
    'real_conditions_exceeding_p99_null': int(np.sum(real_deltas > p99_null_delta)),
    'gate_g5a_verdict': "PASS" if null_false_pos_rate_05 < 0.05 else "FAIL"
}

with open('results/astr/dmso_null_summary.json', 'w') as f:
    json.dump(summary_g5a, f, indent=2)

print("--- Gate G5A Null Calibration Summary ---")
print(json.dumps(summary_g5a, indent=2))
