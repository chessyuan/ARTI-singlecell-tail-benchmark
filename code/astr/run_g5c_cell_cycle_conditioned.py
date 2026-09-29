# code/astr/run_g5c_cell_cycle_conditioned.py
# Gate G5C: Cell-Cycle Conditioned Analysis and Residual Decomposition

import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.stats import chisquare, ttest_rel
import statsmodels.api as sm
from sklearn.neighbors import NearestNeighbors

os.makedirs('tables/astr', exist_ok=True)
os.makedirs('results/astr', exist_ok=True)

print("--- Running Gate G5C: Cell-Cycle Conditioned Analysis ---")
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

# QC filtering
obs['is_qc_excluded'] = (obs['pct_mito'] >= 15.0) | (obs['total_counts'] < 800) | (obs['n_genes'] < 300)
obs['is_qc_pass'] = ~obs['is_dying'] if 'is_dying' in obs.columns else ~obs['is_qc_excluded']
qc_pass_obs = obs[obs['is_qc_pass']].copy()

# Classified states
classified = pd.read_parquet('data_processed/tahoe_mvp_cells_classified.parquet')
qc_pass_obs['cell_state_category'] = classified.loc[qc_pass_obs.index, 'cell_state_category']

# Normalized expression
lib_size = adata.X.sum(axis=1)
if hasattr(lib_size, 'A1'):
    lib_size = lib_size.A1
else:
    lib_size = np.array(lib_size).flatten()
lib_size[lib_size == 0] = 1.0

s_genes = ['MCM2', 'MCM4', 'MCM6', 'PCNA', 'TYMS', 'FEN1', 'RRM1', 'RRM2', 'GINS2']
g2m_genes = ['MKI67', 'TOP2A', 'CDK1', 'CCNB1', 'CCNB2', 'CENPF', 'NUSAP1', 'AURKA', 'BUB1']

var_names = list(adata.var_names)
s_idx = [var_names.index(g) for g in s_genes if g in var_names]
g2m_idx = [var_names.index(g) for g in g2m_genes if g in var_names]

# Compute phase scores
s_sub = adata.X[:, s_idx].toarray()
g2m_sub = adata.X[:, g2m_idx].toarray()

s_norm = np.log1p((s_sub / lib_size[:, None]) * 10000.0).mean(axis=1)
g2m_norm = np.log1p((g2m_sub / lib_size[:, None]) * 10000.0).mean(axis=1)

# Center scores
s_score = s_norm - np.median(s_norm)
g2m_score = g2m_norm - np.median(g2m_norm)

qc_pass_obs['s_score'] = s_score[qc_pass_obs['cell_idx']]
qc_pass_obs['g2m_score'] = g2m_score[qc_pass_obs['cell_idx']]

# Assign phase
phases = []
for s_sc, g2m_sc in zip(qc_pass_obs['s_score'], qc_pass_obs['g2m_score']):
    if s_sc <= 0 and g2m_sc <= 0:
        phases.append('G1')
    elif s_sc > g2m_sc:
        phases.append('S')
    else:
        phases.append('G2M')
qc_pass_obs['assigned_phase'] = phases

# 1. Phase Distribution in AST vs Dominant Response
ast_cells = qc_pass_obs[qc_pass_obs['cell_state_category'] == 'Adaptive_Survival_Tail']
dom_cells = qc_pass_obs[qc_pass_obs['cell_state_category'] == 'Dominant_Response']

ast_phase_counts = ast_cells['assigned_phase'].value_counts()
dom_phase_counts = dom_cells['assigned_phase'].value_counts()
ast_phase_dist = (ast_phase_counts / len(ast_cells)).round(4).to_dict()
dom_phase_dist = (dom_phase_counts / len(dom_cells)).round(4).to_dict()

print(f"AST Phase Distribution: {ast_phase_dist}")
print(f"Dominant Response Phase Distribution: {dom_phase_dist}")

# Chi-square test on phase independence
expected_ast_counts = np.array([dom_phase_dist[p] * len(ast_cells) for p in ['G1', 'S', 'G2M']])
observed_ast_counts = np.array([ast_phase_counts.get(p, 0) for p in ['G1', 'S', 'G2M']])
chi2_stat, chi2_p = chisquare(f_obs=observed_ast_counts, f_exp=expected_ast_counts)

# 2. Residual Regression Analysis: Does AST remain significant after controlling for continuous S & G2M scores?
# Primary model: HOP62 + Adagrasib
hop62_adag = qc_pass_obs[(qc_pass_obs['cell_name'] == 'HOP62') & (qc_pass_obs['drug'] == 'Adagrasib')].copy()
hop62_adag['is_ast'] = (hop62_adag['cell_state_category'] == 'Adaptive_Survival_Tail').astype(int)

# Top 1000 HVGs PCA for HOP62 Adagrasib
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[hop62_adag['cell_idx'], :][:, top_gene_indices].toarray()
sub_lib = sub_X.sum(axis=1, keepdims=True)
sub_lib[sub_lib == 0] = 1.0
sub_norm = np.log1p((sub_X / sub_lib) * 10000.0)

from sklearn.decomposition import PCA
pca = PCA(n_components=5, random_state=42)
pcs = pca.fit_transform(sub_norm)
dists = np.linalg.norm(pcs - np.median(pcs, axis=0, keepdims=True), axis=1)
hop62_adag['dispersion_distance'] = dists

# OLS Regression: dispersion_distance ~ is_ast + s_score + g2m_score
X_reg = sm.add_constant(hop62_adag[['is_ast', 's_score', 'g2m_score']])
y_reg = hop62_adag['dispersion_distance']
ols_model = sm.OLS(y_reg, X_reg).fit()

ast_coef = float(ols_model.params['is_ast'])
ast_pval = float(ols_model.pvalues['is_ast'])
r2_val = float(ols_model.rsquared)

print(f"OLS Regression in HOP62 Adagrasib (dispersion ~ is_ast + s_score + g2m_score):")
print(f"  AST coefficient: {ast_coef:.4f}, p-value: {ast_pval:.2e}, R2: {r2_val:.4f}")

# 3. Cell-Cycle Matched Nearest-Neighbor Comparison
# Match each AST cell to 1 nearest Dominant cell based on (s_score, g2m_score)
ast_sub = hop62_adag[hop62_adag['is_ast'] == 1]
dom_sub = hop62_adag[hop62_adag['is_ast'] == 0]

nn = NearestNeighbors(n_neighbors=1)
nn.fit(dom_sub[['s_score', 'g2m_score']])
distances, indices = nn.kneighbors(ast_sub[['s_score', 'g2m_score']])

matched_dom_dists = dom_sub['dispersion_distance'].iloc[indices.flatten()].values
ast_dists = ast_sub['dispersion_distance'].values

ttest_stat, ttest_pval = ttest_rel(ast_dists, matched_dom_dists)
mean_ast_dist = float(np.mean(ast_dists))
mean_matched_dom_dist = float(np.mean(matched_dom_dists))

print(f"Cell-Cycle Matched Comparison:")
print(f"  Mean AST distance: {mean_ast_dist:.4f}")
print(f"  Mean Matched Dominant distance: {mean_matched_dom_dist:.4f}")
print(f"  Paired t-test: t = {ttest_stat:.4f}, p = {ttest_pval:.4f}")

summary_g5c = {
    'ast_cells_evaluated': len(ast_cells),
    'dominant_cells_evaluated': len(dom_cells),
    'ast_phase_distribution': ast_phase_dist,
    'dominant_phase_distribution': dom_phase_dist,
    'phase_chi2_statistic': round(float(chi2_stat), 4),
    'phase_chi2_pvalue': float(chi2_p),
    'hop62_adagrasib_ols_ast_coef': round(ast_coef, 4),
    'hop62_adagrasib_ols_ast_pvalue': float(ast_pval),
    'hop62_adagrasib_ols_r2': round(r2_val, 4),
    'matched_ast_mean_distance': round(mean_ast_dist, 4),
    'matched_dominant_mean_distance': round(mean_matched_dom_dist, 4),
    'paired_ttest_pvalue': float(ttest_pval),
    'cell_cycle_independence_verdict': "PASS" if ast_pval < 0.01 and ttest_pval < 0.05 else "FAIL"
}

with open('results/astr/cell_cycle_conditioning_summary.json', 'w') as f:
    json.dump(summary_g5c, f, indent=2)

# Save condition-level phase distribution table
cond_phase_list = []
for (cname, drug), grp in qc_pass_obs[qc_pass_obs['drug'] != 'DMSO_TF'].groupby(['cell_name', 'drug'], observed=True):
    if len(grp) >= 25:
        c_ast = grp[grp['cell_state_category'] == 'Adaptive_Survival_Tail']
        cond_phase_list.append({
            'cell_name': cname,
            'drug': drug,
            'n_cells': len(grp),
            'n_ast': len(c_ast),
            'g1_fraction': round(float((grp['assigned_phase'] == 'G1').mean()), 3),
            's_fraction': round(float((grp['assigned_phase'] == 'S').mean()), 3),
            'g2m_fraction': round(float((grp['assigned_phase'] == 'G2M').mean()), 3),
            'ast_g1_fraction': round(float((c_ast['assigned_phase'] == 'G1').mean()), 3) if len(c_ast) > 0 else np.nan,
            'ast_s_fraction': round(float((c_ast['assigned_phase'] == 'S').mean()), 3) if len(c_ast) > 0 else np.nan,
            'ast_g2m_fraction': round(float((c_ast['assigned_phase'] == 'G2M').mean()), 3) if len(c_ast) > 0 else np.nan
        })

pd.DataFrame(cond_phase_list).to_csv('tables/astr/ast_cell_cycle_conditioned.csv', index=False)
print("Saved tables/astr/ast_cell_cycle_conditioned.csv successfully!")
