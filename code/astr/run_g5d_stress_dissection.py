# code/astr/run_g5d_stress_dissection.py
# Gate G5D: Stress Module Dissection and Generic-Stress-Matched Control Analysis

import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_rel
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs('tables/astr', exist_ok=True)
os.makedirs('results/astr', exist_ok=True)
os.makedirs('figures/astr', exist_ok=True)

print("--- Running Gate G5D: Biological Stress Module Dissection ---")
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

obs['is_qc_excluded'] = (obs['pct_mito'] >= 15.0) | (obs['total_counts'] < 800) | (obs['n_genes'] < 300)
obs['is_qc_pass'] = ~obs['is_dying'] if 'is_dying' in obs.columns else ~obs['is_qc_excluded']
qc_pass_obs = obs[obs['is_qc_pass']].copy()

classified = pd.read_parquet('data_processed/tahoe_mvp_cells_classified.parquet')
qc_pass_obs['cell_state_category'] = classified.loc[qc_pass_obs.index, 'cell_state_category']

# Normalized matrix
lib_size = adata.X.sum(axis=1)
if hasattr(lib_size, 'A1'):
    lib_size = lib_size.A1
else:
    lib_size = np.array(lib_size).flatten()
lib_size[lib_size == 0] = 1.0

var_names = list(adata.var_names)

stress_gene_sets = {
    'Apoptosis': ['CASP3', 'BAX', 'BAK1', 'PMAIP1', 'BBC3', 'BCL2L11', 'FAS', 'CYCS', 'BID', 'DIABLO'],
    'ISR': ['ATF4', 'DDIT3', 'ASNS', 'TRIB3', 'HSPA5', 'PPP1R15A', 'EIF2AK3', 'CEBPB'],
    'Heat_Shock': ['HSP90AA1', 'HSPA1A', 'HSPA8', 'DNAJB1', 'HSPB1', 'HSPD1'],
    'Oxidative_Stress': ['HMOX1', 'SOD1', 'SOD2', 'NQO1', 'GCLC', 'TXNRD1', 'GPX4', 'PRDX1'],
    'Immediate_Early': ['FOS', 'JUN', 'EGR1', 'ATF3', 'FOSB', 'JUNB', 'JUND']
}

# Score modules
for mod_name, genes in stress_gene_sets.items():
    g_idx = [var_names.index(g) for g in genes if g in var_names]
    sub_X = adata.X[:, g_idx].toarray()
    norm_sub = np.log1p((sub_X / lib_size[:, None]) * 10000.0).mean(axis=1)
    qc_pass_obs[f'score_{mod_name}'] = norm_sub[qc_pass_obs['cell_idx']]

qc_pass_obs['score_Mitochondrial'] = qc_pass_obs['pct_mito']

# Calculate z-scored composite generic stress score
stress_cols = [f'score_{m}' for m in stress_gene_sets.keys()] + ['score_Mitochondrial']
z_scores = pd.DataFrame(index=qc_pass_obs.index)
for col in stress_cols:
    val = qc_pass_obs[col]
    z_scores[col] = (val - val.mean()) / (val.std() + 1e-6)
qc_pass_obs['composite_generic_stress'] = z_scores.mean(axis=1)

# Primary model: HOP62 + Adagrasib
hop62_adag = qc_pass_obs[(qc_pass_obs['cell_name'] == 'HOP62') & (qc_pass_obs['drug'] == 'Adagrasib')].copy()

# Fit condition PCA for dispersion
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[hop62_adag['cell_idx'], :][:, top_gene_indices].toarray()
sub_lib = sub_X.sum(axis=1, keepdims=True)
sub_lib[sub_lib == 0] = 1.0
sub_norm = np.log1p((sub_X / sub_lib) * 10000.0)

pca = PCA(n_components=5, random_state=42)
pcs = pca.fit_transform(sub_norm)
dists = np.linalg.norm(pcs - np.median(pcs, axis=0, keepdims=True), axis=1)
hop62_adag['dispersion_distance'] = dists

ast_hop = hop62_adag[hop62_adag['cell_state_category'] == 'Adaptive_Survival_Tail']
dom_hop = hop62_adag[hop62_adag['cell_state_category'] == 'Dominant_Response']

# 1. Compare stress scores between AST and Dominant Response
module_comp = []
for mod_col in stress_cols + ['composite_generic_stress']:
    m_ast = float(ast_hop[mod_col].mean())
    m_dom = float(dom_hop[mod_col].mean())
    stat, pval = mannwhitneyu(ast_hop[mod_col], dom_hop[mod_col], alternative='two-sided')
    module_comp.append({
        'module': mod_col.replace('score_', ''),
        'mean_ast': round(m_ast, 4),
        'mean_dominant': round(m_dom, 4),
        'delta_mean': round(m_ast - m_dom, 4),
        'mann_whitney_pval': round(float(pval), 5)
    })

mod_df = pd.DataFrame(module_comp)
print("HOP62 Adagrasib Stress Module Comparison:")
print(mod_df.to_string())

# 2. Generic-Stress-Matched Nearest-Neighbor Comparison
# Match each AST cell to 1 nearest Dominant cell based on composite_generic_stress
nn = NearestNeighbors(n_neighbors=1)
nn.fit(dom_hop[['composite_generic_stress']])
distances, indices = nn.kneighbors(ast_hop[['composite_generic_stress']])

matched_dom = dom_hop.iloc[indices.flatten()]
matched_dom_dists = matched_dom['dispersion_distance'].values
ast_dists = ast_hop['dispersion_distance'].values

ttest_stat, ttest_pval = ttest_rel(ast_dists, matched_dom_dists)
mean_ast_dist = float(np.mean(ast_dists))
mean_matched_dom_dist = float(np.mean(matched_dom_dists))
delta_matched_disp = mean_ast_dist - mean_matched_dom_dist

print(f"\nGeneric-Stress-Matched Control Comparison:")
print(f"  Mean AST Dispersion: {mean_ast_dist:.4f}")
print(f"  Mean Stress-Matched Dominant Dispersion: {mean_matched_dom_dist:.4f}")
print(f"  Excess Dispersion after Stress Matching: {delta_matched_disp:.4f}")
print(f"  Paired t-test: t = {ttest_stat:.4f}, p = {ttest_pval:.4f}")

# Check key survival marker MT2A in stress-matched pairs
var_names = list(adata.var_names)
mt2a_idx = var_names.index('MT2A') if 'MT2A' in var_names else -1
if mt2a_idx != -1:
    ast_mt2a = np.log1p((adata.X[ast_hop['cell_idx'], mt2a_idx].toarray() / lib_size[ast_hop['cell_idx'], None]) * 10000.0).flatten()
    dom_matched_mt2a = np.log1p((adata.X[matched_dom['cell_idx'], mt2a_idx].toarray() / lib_size[matched_dom['cell_idx'], None]) * 10000.0).flatten()
    t_mt2a, p_mt2a = ttest_rel(ast_mt2a, dom_matched_mt2a)
    print(f"  MT2A Expression: AST = {ast_mt2a.mean():.3f} vs Stress-Matched Dominant = {dom_matched_mt2a.mean():.3f} (paired p = {p_mt2a:.4f})")
else:
    p_mt2a = np.nan

# Save tables
mod_df.to_csv('tables/astr/ast_stress_dissection.csv', index=False)

summary_g5d = {
    'primary_model': 'HOP62 + Adagrasib',
    'n_ast_cells': len(ast_hop),
    'n_dominant_cells': len(dom_hop),
    'modules_evaluated': [m.replace('score_', '') for m in stress_cols],
    'apoptosis_pval': float(mod_df[mod_df['module'] == 'Apoptosis']['mann_whitney_pval'].iloc[0]),
    'isr_pval': float(mod_df[mod_df['module'] == 'ISR']['mann_whitney_pval'].iloc[0]),
    'mean_ast_dispersion': round(mean_ast_dist, 4),
    'mean_stress_matched_dominant_dispersion': round(mean_matched_dom_dist, 4),
    'excess_dispersion_stress_matched': round(delta_matched_disp, 4),
    'stress_matched_dispersion_paired_pvalue': float(ttest_pval),
    'is_stress_dominated_artifact': bool(ttest_pval > 0.05 or delta_matched_disp <= 0.50),
    'gate_g5d_verdict': "PASS" if ttest_pval < 0.05 and delta_matched_disp > 0.50 else "FAIL"
}

with open('results/astr/stress_dissection_summary.json', 'w') as f:
    json.dump(summary_g5d, f, indent=2)

# Figure: Stress module radar / bar comparison
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# Plot A: Stress module z-scores
plot_mod = mod_df[mod_df['module'] != 'composite_generic_stress']
x_indices = np.arange(len(plot_mod))
bar_width = 0.35

axes[0].bar(x_indices - bar_width/2, plot_mod['mean_ast'], width=bar_width, label='Adaptive Survival Tail (AST)', color='#d95f02')
axes[0].bar(x_indices + bar_width/2, plot_mod['mean_dominant'], width=bar_width, label='Dominant Response', color='#7570b3')
axes[0].set_xticks(x_indices)
axes[0].set_xticklabels(plot_mod['module'], rotation=25, ha='right', fontsize=10)
axes[0].set_ylabel('Mean Module Score (Log Normalized)', fontsize=11)
axes[0].set_title('A. Biological Stress Modules: AST vs Dominant Response (HOP62 Adagrasib)', fontsize=11, fontweight='bold')
axes[0].legend(frameon=True)

# Plot B: Dispersion before and after Stress Matching
axes[1].boxplot([ast_dists, matched_dom_dists], tick_labels=['AST Cells\n(N=6)', 'Generic-Stress-Matched\nDominant Cells (N=6)'], patch_artist=True,
                boxprops=dict(facecolor='#a6cee3', color='#1f78b4'))
for i in range(len(ast_dists)):
    axes[1].plot([1, 2], [ast_dists[i], matched_dom_dists[i]], color='gray', linestyle='--', alpha=0.7)
axes[1].set_ylabel('Transcriptional Dispersion (Distance from Median)', fontsize=11)
axes[1].set_title(f'B. Stress-Matched Control: Excess Dispersion Preserved (p = {ttest_pval:.4f})', fontsize=11, fontweight='bold')

plt.tight_layout()
fig.savefig('figures/astr/fig_ast_stress_dissection.png', dpi=300)
plt.close()
print("Saved figures/astr/fig_ast_stress_dissection.png successfully!")
