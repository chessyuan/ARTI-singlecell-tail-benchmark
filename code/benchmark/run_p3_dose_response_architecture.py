import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, ttest_rel
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs('tables/benchmark', exist_ok=True)
os.makedirs('results/benchmark', exist_ok=True)
os.makedirs('figures/benchmark', exist_ok=True)

print("--- Running P3: Dose-Response Architecture Analysis ---")

adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

obs['is_qc_excluded'] = (obs['pct_mito'] >= 15.0) | (obs['total_counts'] < 800) | (obs['n_genes'] < 300)
obs['is_qc_pass'] = ~obs['is_dying'] if 'is_dying' in obs.columns else ~obs['is_qc_excluded']
qc_pass_obs = obs[obs['is_qc_pass']].copy()

# 1000 HVGs
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[:, top_gene_indices].toarray()
lib_size = sub_X.sum(axis=1, keepdims=True)
lib_size[lib_size == 0] = 1.0
norm_X = np.log1p((sub_X / lib_size) * 10000.0)

from sklearn.decomposition import PCA
pca = PCA(n_components=5, random_state=42)
pca_X = pca.fit_transform(norm_X)

# Precompute DMSO data per cell line and plate
dmso_profiles = {}
for (cname, plate), grp in qc_pass_obs[qc_pass_obs['drug'] == 'DMSO_TF'].groupby(['cell_name', 'plate'], observed=True):
    idx = grp['cell_idx'].values
    dmso_profiles[(cname, plate)] = {
        'norm_mean': norm_X[idx].mean(axis=0),
        'pca_X': pca_X[idx],
        'pca_mean': pca_X[idx].mean(axis=0),
        'pca_dispersion': float(np.mean(np.linalg.norm(pca_X[idx] - pca_X[idx].mean(axis=0), axis=1))),
        'count': len(idx)
    }

# Find cell_name + drug with both 0.05 uM and 5.0 uM across plates
pvt_counts = qc_pass_obs[qc_pass_obs['drug'] != 'DMSO_TF'].groupby(['cell_name', 'drug', 'dose'], observed=True).size().unstack(fill_value=0)
valid_pairs = pvt_counts[(pvt_counts['0.05 uM'] >= 10) & (pvt_counts['5.0 uM'] >= 10)].index

dose_records = []

for cname, drug in valid_pairs:
    # Filter cells for low and high dose
    low_cells = qc_pass_obs[(qc_pass_obs['cell_name'] == cname) & (qc_pass_obs['drug'] == drug) & (qc_pass_obs['dose'] == '0.05 uM')]
    high_cells = qc_pass_obs[(qc_pass_obs['cell_name'] == cname) & (qc_pass_obs['drug'] == drug) & (qc_pass_obs['dose'] == '5.0 uM')]
    
    # Use dominant plate for matching DMSO
    low_plate = low_cells['plate'].value_counts().index[0]
    high_plate = high_cells['plate'].value_counts().index[0]
    
    dmso_low = dmso_profiles.get((cname, low_plate))
    dmso_high = dmso_profiles.get((cname, high_plate))
    
    if dmso_low is None or dmso_high is None:
        continue
        
    low_idx = low_cells['cell_idx'].values
    high_idx = high_cells['cell_idx'].values
    
    # Mean response (AR1 pseudobulk L2)
    low_norm_mean = norm_X[low_idx].mean(axis=0)
    high_norm_mean = norm_X[high_idx].mean(axis=0)
    ar1_low = float(np.linalg.norm(low_norm_mean - dmso_low['norm_mean']))
    ar1_high = float(np.linalg.norm(high_norm_mean - dmso_high['norm_mean']))
    
    # ARTI (dispersion relative to DMSO / 2.0)
    low_pca = pca_X[low_idx]
    high_pca = pca_X[high_idx]
    
    low_disp = float(np.mean(np.linalg.norm(low_pca - low_pca.mean(axis=0), axis=1)))
    high_disp = float(np.mean(np.linalg.norm(high_pca - high_pca.mean(axis=0), axis=1)))
    
    arti_low = float((low_disp - dmso_low['pca_dispersion']) / 2.0)
    arti_high = float((high_disp - dmso_high['pca_dispersion']) / 2.0)
    
    delta_mean = ar1_high - ar1_low
    delta_arti = arti_high - arti_low
    
    # Architecture classification
    if delta_mean > 0 and delta_arti > 0:
        arch = 'Monotonic_Co_Expansion'
    elif delta_mean > 0 and delta_arti < 0:
        arch = 'Adaptive_Heterogeneity_Window'
    elif delta_mean <= 0 and delta_arti > 0:
        arch = 'Tail_Shift_Without_Bulk_Effect'
    else:
        arch = 'Sub_Threshold_or_Tolerant'
        
    dose_records.append({
        'cell_name': cname,
        'drug': drug,
        'n_cells_low': len(low_cells),
        'n_cells_high': len(high_cells),
        'ar1_mean_response_low': ar1_low,
        'ar1_mean_response_high': ar1_high,
        'delta_mean_response': delta_mean,
        'arti_low': arti_low,
        'arti_high': arti_high,
        'delta_arti': delta_arti,
        'architecture_category': arch
    })

dose_df = pd.DataFrame(dose_records)
dose_df.to_csv('tables/benchmark/dose_response_architecture.csv', index=False)
print(f"Evaluated dose architecture across {len(dose_df)} drug-cell line pairs.")

# Summary statistics
arch_counts = dose_df['architecture_category'].value_counts()
print("\nDose-Response Architecture Breakdown:")
for cat, cnt in arch_counts.items():
    print(f"  {cat:32s}: {cnt:3d} ({cnt/len(dose_df)*100:.1f}%)")

# Paired test on Mean Response and ARTI
t_mean, p_mean = ttest_rel(dose_df['ar1_mean_response_high'], dose_df['ar1_mean_response_low'])
t_arti, p_arti = ttest_rel(dose_df['arti_high'], dose_df['arti_low'])

print(f"\nPaired Test across 100-Fold Dose Escalation (0.05 uM -> 5.0 uM):")
print(f"  Mean Response (AR1): Low = {dose_df['ar1_mean_response_low'].mean():.3f}, High = {dose_df['ar1_mean_response_high'].mean():.3f} (t = {t_mean:+.2f}, p = {p_mean:.2e})")
print(f"  Tail Index (ARTI):   Low = {dose_df['arti_low'].mean():.3f}, High = {dose_df['arti_high'].mean():.3f} (t = {t_arti:+.2f}, p = {p_arti:.2e})")

# Correlation between Delta Mean and Delta ARTI
pr_delta, p_delta = pearsonr(dose_df['delta_mean_response'], dose_df['delta_arti'])
print(f"  Correlation between Delta Mean and Delta ARTI: r = {pr_delta:+.4f} (p = {p_delta:.2e})")

# Permutation test on dose order
np.random.seed(42)
n_perms = 1000
obs_diff_arti = dose_df['delta_arti'].mean()
perm_diffs = []
for _ in range(n_perms):
    # Flip signs randomly
    signs = np.random.choice([-1, 1], size=len(dose_df))
    perm_diffs.append(np.mean(dose_df['delta_arti'] * signs))

emp_p_dose = float(np.mean(np.abs(perm_diffs) >= np.abs(obs_diff_arti)))
print(f"  Dose Order Permutation Test: Observed Diff = {obs_diff_arti:+.4f}, Empirical p = {emp_p_dose:.4f}")

dose_summary = {
    'pairs_evaluated': len(dose_df),
    'architecture_counts': arch_counts.to_dict(),
    'architecture_fractions': (arch_counts / len(dose_df)).to_dict(),
    'paired_t_tests': {
        'mean_response_low': float(dose_df['ar1_mean_response_low'].mean()),
        'mean_response_high': float(dose_df['ar1_mean_response_high'].mean()),
        'mean_response_p': float(p_mean),
        'arti_low': float(dose_df['arti_low'].mean()),
        'arti_high': float(dose_df['arti_high'].mean()),
        'arti_p': float(p_arti)
    },
    'correlation_delta_mean_vs_delta_arti': {
        'pearson_r': float(pr_delta),
        'p_value': float(p_delta)
    },
    'dose_order_permutation_p': float(emp_p_dose)
}

with open('results/benchmark/dose_response_summary.json', 'w') as f:
    json.dump(dose_summary, f, indent=2)

# Figure: Multi-panel dose response visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Panel A: Mean Response vs ARTI slope
palette = {
    'Monotonic_Co_Expansion': '#2166ac',
    'Adaptive_Heterogeneity_Window': '#b2182b',
    'Tail_Shift_Without_Bulk_Effect': '#1b7837',
    'Sub_Threshold_or_Tolerant': '#999999'
}

sns.scatterplot(data=dose_df, x='delta_mean_response', y='delta_arti',
                hue='architecture_category', palette=palette, s=55, alpha=0.8, ax=ax1)
ax1.axhline(0, color='gray', linestyle='--', alpha=0.6)
ax1.axvline(0, color='gray', linestyle='--', alpha=0.6)
ax1.set_title('A. Dose-Induced Changes: Delta Mean Response vs. Delta ARTI\n(100-Fold Concentration Escalation: 0.05 uM to 5.0 uM)', fontsize=11)
ax1.set_xlabel('Delta Mean Response (High - Low AR1 L2 Distance)', fontsize=10)
ax1.set_ylabel('Delta ARTI (High - Low Tail Index)', fontsize=10)
ax1.grid(True, linestyle='--', alpha=0.4)
ax1.legend(title='Architecture', loc='upper left', fontsize=8)

# Panel B: Boxplots comparing Low vs High Dose
melted = pd.melt(dose_df[['ar1_mean_response_low', 'ar1_mean_response_high', 'arti_low', 'arti_high']],
                 var_name='condition', value_name='score')
melted['metric'] = melted['condition'].apply(lambda x: 'Mean Response (AR1)' if 'ar1' in x else 'Tail Index (ARTI)')
melted['dose'] = melted['condition'].apply(lambda x: '0.05 uM (Low)' if 'low' in x else '5.0 uM (High)')

sns.boxplot(data=melted, x='metric', y='score', hue='dose', palette=['#67a9cf', '#ef8a62'], ax=ax2, width=0.5)
ax2.set_title(f'B. Population Level Comparison (N = {len(dose_df)} Pairs)\nMean Shift increases (p = {p_mean:.1e}) while ARTI exhibits non-uniform divergence', fontsize=11)
ax2.set_xlabel('')
ax2.set_ylabel('Standardized Score', fontsize=10)
ax2.grid(True, linestyle='--', alpha=0.4)
ax2.legend(title='Dose Level', loc='upper right', fontsize=9)

plt.tight_layout()
plt.savefig('figures/benchmark/fig_dose_response_architecture.png', dpi=300)
plt.close()
print("Saved figure to figures/benchmark/fig_dose_response_architecture.png")
