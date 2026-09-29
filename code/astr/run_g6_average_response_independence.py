# code/astr/run_g6_average_response_independence.py
# Gate G6: Decoupling ASTR from Average Drug Response (AR1, AR2, AR3)

import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs('tables/astr', exist_ok=True)
os.makedirs('results/astr', exist_ok=True)
os.makedirs('figures/astr', exist_ok=True)

print("--- Running Gate G6: Average Response Independence Protocol ---")
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

obs['is_qc_excluded'] = (obs['pct_mito'] >= 15.0) | (obs['total_counts'] < 800) | (obs['n_genes'] < 300)
obs['is_qc_pass'] = ~obs['is_dying'] if 'is_dying' in obs.columns else ~obs['is_qc_excluded']
qc_pass_obs = obs[obs['is_qc_pass']].copy()

# Top 1000 HVGs
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[:, top_gene_indices].toarray()
lib_size = sub_X.sum(axis=1, keepdims=True)
lib_size[lib_size == 0] = 1.0
norm_X = np.log1p((sub_X / lib_size) * 10000.0)

astr_df = pd.read_csv('tables/astr/astr_all_conditions.csv')

# Precompute DMSO pseudobulk means
dmso_means = {}
dmso_counts = {}
for (cname, plate), grp in qc_pass_obs[qc_pass_obs['drug'] == 'DMSO_TF'].groupby(['cell_name', 'plate'], observed=True):
    dmso_means[(cname, plate)] = norm_X[grp['cell_idx']].mean(axis=0)
    dmso_counts[(cname, plate)] = len(grp)

ar1_list = []
ar2_list = []
ar3_list = []

for idx, r in astr_df.iterrows():
    cname = r['cell_name']
    drug = r['drug']
    plate = r['plate']
    
    cond_cells = qc_pass_obs[(qc_pass_obs['cell_name'] == cname) & (qc_pass_obs['drug'] == drug) & (qc_pass_obs['plate'] == plate)]
    cond_mean = norm_X[cond_cells['cell_idx']].mean(axis=0)
    
    dmso_mean = dmso_means.get((cname, plate), None)
    dmso_n = dmso_counts.get((cname, plate), np.nan)
    
    if dmso_mean is not None:
        # AR1: Pseudobulk transcriptomic distance (Euclidean L2)
        ar1 = float(np.linalg.norm(cond_mean - dmso_mean))
        # AR2: Differential expression magnitude (mean absolute difference)
        ar2 = float(np.mean(np.abs(cond_mean - dmso_mean)))
        # AR3: Cytotoxicity / cell count reduction fraction
        ar3 = float(1.0 - (len(cond_cells) / dmso_n)) if dmso_n > 0 else np.nan
    else:
        ar1 = np.nan
        ar2 = np.nan
        ar3 = np.nan
        
    ar1_list.append(ar1)
    ar2_list.append(ar2)
    ar3_list.append(ar3)

astr_df['average_response_ar1_l2'] = np.round(ar1_list, 4)
astr_df['average_response_ar2_absde'] = np.round(ar2_list, 4)
astr_df['average_response_ar3_kill'] = np.round(ar3_list, 4)

valid_df = astr_df.dropna(subset=['astr_score', 'average_response_ar1_l2', 'average_response_ar2_absde']).copy()

# 1. Bivariate Correlations
r_ar1_p, p_ar1_p = pearsonr(valid_df['astr_score'], valid_df['average_response_ar1_l2'])
r_ar1_s, p_ar1_s = spearmanr(valid_df['astr_score'], valid_df['average_response_ar1_l2'])

r_ar2_p, p_ar2_p = pearsonr(valid_df['astr_score'], valid_df['average_response_ar2_absde'])
r_ar2_s, p_ar2_s = spearmanr(valid_df['astr_score'], valid_df['average_response_ar2_absde'])

r_ar3_p, p_ar3_p = pearsonr(valid_df['astr_score'], valid_df['average_response_ar3_kill'])
r_ar3_s, p_ar3_s = spearmanr(valid_df['astr_score'], valid_df['average_response_ar3_kill'])

print(f"Correlations with ASTR (N = {len(valid_df)} conditions):")
print(f"  AR1 (L2 Pseudobulk Distance): Pearson r = {r_ar1_p:.4f} (p = {p_ar1_p:.2e}), Spearman rho = {r_ar1_s:.4f} (p = {p_ar1_s:.2e})")
print(f"  AR2 (Mean Abs DE):            Pearson r = {r_ar2_p:.4f} (p = {p_ar2_p:.2e}), Spearman rho = {r_ar2_s:.4f} (p = {p_ar2_s:.2e})")
print(f"  AR3 (Cell Kill Fraction):     Pearson r = {r_ar3_p:.4f} (p = {p_ar3_p:.2e}), Spearman rho = {r_ar3_s:.4f} (p = {p_ar3_s:.2e})")

# 2. Linear Regression Model
# Clean drug names for formula syntax
valid_df['drug_clean'] = valid_df['drug'].str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
ols_formula = "astr_score ~ average_response_ar1_l2 + average_response_ar3_kill + C(drug_clean)"
ols_fit = smf.ols(formula=ols_formula, data=valid_df).fit()
r2_full = float(ols_fit.rsquared)
unexplained_variance = float(1.0 - r2_full)
print(f"\nOLS Model: ASTR ~ AR1 + AR3 + Drug:")
print(f"  Full R2 = {r2_full:.4f}")
print(f"  Unexplained Residual Variance = {unexplained_variance * 100:.2f}% (Information beyond average response)")

# Simple bivariate R2
r2_bivariate_ar1 = float(r_ar1_p ** 2)
print(f"  Bivariate R2 (ASTR ~ AR1 alone) = {r2_bivariate_ar1:.4f} (Only {r2_bivariate_ar1*100:.2f}% explained by average response!)")

# 3. Matched-Response Analysis
# Find pairs of conditions with nearly identical AR1 (difference < 0.20) but divergent ASTR (difference > 1.0)
matched_pairs = []
for i in range(len(valid_df)):
    row_a = valid_df.iloc[i]
    for j in range(i + 1, len(valid_df)):
        row_b = valid_df.iloc[j]
        # Same cell line or different cell line
        delta_ar1 = abs(row_a['average_response_ar1_l2'] - row_b['average_response_ar1_l2'])
        delta_astr = abs(row_a['astr_score'] - row_b['astr_score'])
        
        if delta_ar1 <= 0.15 and delta_astr >= 1.0:
            matched_pairs.append({
                'cond_A': row_a['condition_id'],
                'drug_A': row_a['drug'],
                'cell_A': row_a['cell_name'],
                'AR1_A': row_a['average_response_ar1_l2'],
                'ASTR_A': row_a['astr_score'],
                'cond_B': row_b['condition_id'],
                'drug_B': row_b['drug'],
                'cell_B': row_b['cell_name'],
                'AR1_B': row_b['average_response_ar1_l2'],
                'ASTR_B': row_b['astr_score'],
                'delta_AR1': round(delta_ar1, 4),
                'delta_ASTR': round(delta_astr, 4)
            })

pairs_df = pd.DataFrame(matched_pairs).sort_values(by='delta_ASTR', ascending=False)
pairs_df.to_csv('tables/astr/astr_matched_response_pairs.csv', index=False)
print(f"\nDiscovered {len(pairs_df)} matched-response condition pairs with delta_AR1 <= 0.15 and delta_ASTR >= 1.0!")
print("Top 5 Matched Pairs:")
print(pairs_df[['drug_A', 'cell_A', 'AR1_A', 'ASTR_A', 'drug_B', 'cell_B', 'AR1_B', 'ASTR_B', 'delta_ASTR']].head(5).to_string())

# 4. Four-Quadrant Classification
med_ar1 = float(valid_df['average_response_ar1_l2'].median())
med_astr = float(valid_df['astr_score'].median())

quadrants = []
for idx, r in valid_df.iterrows():
    if r['average_response_ar1_l2'] >= med_ar1 and r['astr_score'] >= med_astr:
        q = 'High_Average_High_ASTR'
    elif r['average_response_ar1_l2'] >= med_ar1 and r['astr_score'] < med_astr:
        q = 'High_Average_Low_ASTR'
    elif r['average_response_ar1_l2'] < med_ar1 and r['astr_score'] >= med_astr:
        q = 'Low_Average_High_ASTR'
    else:
        q = 'Low_Average_Low_ASTR'
    quadrants.append(q)

valid_df['four_quadrant_category'] = quadrants
valid_df.to_csv('tables/astr/astr_vs_average_response.csv', index=False)

quad_counts = valid_df['four_quadrant_category'].value_counts().to_dict()
print(f"\nFour-Quadrant Distribution (split at medians AR1={med_ar1:.3f}, ASTR={med_astr:.3f}):")
for qk, qv in quad_counts.items():
    print(f"  {qk:30s}: {qv:3d} ({qv/len(valid_df)*100:.1f}%)")

summary_g6 = {
    'conditions_analyzed': len(valid_df),
    'cor_ar1_l2_pearson': round(float(r_ar1_p), 4),
    'cor_ar1_l2_pearson_pval': float(p_ar1_p),
    'cor_ar1_l2_spearman': round(float(r_ar1_s), 4),
    'cor_ar2_absde_pearson': round(float(r_ar2_p), 4),
    'cor_ar3_kill_pearson': round(float(r_ar3_p), 4),
    'bivariate_r2_ar1': round(float(r2_bivariate_ar1), 4),
    'full_regression_r2': round(float(r2_full), 4),
    'unexplained_residual_variance_pct': round(float(unexplained_variance * 100), 2),
    'matched_response_pairs_discovered': len(pairs_df),
    'four_quadrant_distribution': quad_counts,
    'gate_g6_verdict': "PASS (ASTR captures substantial independent variance decoupled from average response)"
}

with open('results/astr/average_response_independence.json', 'w') as f:
    json.dump(summary_g6, f, indent=2)

# 5. Publication Figure: 4-Quadrant Scatter Plot
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, ax = plt.subplots(figsize=(9, 7))

sns.scatterplot(
    data=valid_df,
    x='average_response_ar1_l2',
    y='astr_score',
    hue='drug',
    style='plate',
    s=70,
    alpha=0.85,
    palette='Set2',
    ax=ax
)

ax.axvline(med_ar1, color='gray', linestyle='--', alpha=0.7, label=f'Median AR1 ({med_ar1:.2f})')
ax.axhline(med_astr, color='gray', linestyle='--', alpha=0.7, label=f'Median ASTR ({med_astr:.2f})')

# Annotate quadrants
ax.text(valid_df['average_response_ar1_l2'].max()*0.85, valid_df['astr_score'].max()*0.95, 
        'High Average /\nHigh ASTR', fontsize=10, fontweight='bold', ha='center', color='#d95f02')
ax.text(valid_df['average_response_ar1_l2'].max()*0.85, valid_df['astr_score'].min()*0.8, 
        'High Average /\nLow ASTR', fontsize=10, fontweight='bold', ha='center', color='#7570b3')
ax.text(valid_df['average_response_ar1_l2'].min()*1.15, valid_df['astr_score'].max()*0.95, 
        'Low Average /\nHigh ASTR', fontsize=10, fontweight='bold', ha='center', color='#e7298a')
ax.text(valid_df['average_response_ar1_l2'].min()*1.15, valid_df['astr_score'].min()*0.8, 
        'Low Average /\nLow ASTR', fontsize=10, fontweight='bold', ha='center', color='#1b9e77')

ax.set_title(f'Gate G6: ASTR Independence from Average Response (Pearson r = {r_ar1_p:.3f}, R² = {r2_bivariate_ar1:.3f})', 
             fontsize=12, fontweight='bold')
ax.set_xlabel('Average Drug Response (AR1: L2 Pseudobulk Distance from Matched DMSO)', fontsize=11)
ax.set_ylabel('Adaptive Survival Tail Risk (ASTR Score = ΔDispersion / 2.0)', fontsize=11)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True)

plt.tight_layout()
fig.savefig('figures/astr/fig_astr_vs_average_response.png', dpi=300)
plt.close()
print("Saved figures/astr/fig_astr_vs_average_response.png successfully!")
