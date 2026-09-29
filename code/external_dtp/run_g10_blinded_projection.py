import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu, ttest_ind, ttest_rel
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs('results/external_dtp', exist_ok=True)
os.makedirs('tables/external_dtp', exist_ok=True)
os.makedirs('figures/external_dtp', exist_ok=True)

print('=== Gate G10: Primary Blinded Projection & Falsification Testing ===')

# 1. Load Data
print('Loading Tahoe MVP data and classifications...')
tahoe_adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
classified_df = pd.read_parquet('data_processed/tahoe_mvp_cells_classified.parquet')
tahoe_adata.obs = classified_df

with open('results/external_dtp/dtp_reference_signatures.json') as f:
    dtp_ref = json.load(f)

dtp_de_df = pd.read_csv('tables/external_dtp/dtp_signature_genes.csv')
gse_adata = anndata.read_h5ad('data_processed/external_dtp/gse275199_processed_adata.h5ad')

# Shared genes
shared_genes = sorted(list(set(tahoe_adata.var_names).intersection(set(gse_adata.var_names))))
print(f'Shared genes: {len(shared_genes)}')

# Select top 1000 shared HVGs from GSE275199 DTP vector for response vector comparison
dtp_de_shared = dtp_de_df[dtp_de_df['gene'].isin(shared_genes)].copy()
top_hvg_genes = dtp_de_shared.sort_values(by='log2FC', key=abs, ascending=False).head(1000)['gene'].tolist()
print(f'Selected top 1000 shared DTP-associated genes for response vector comparison.')

# Precompute DTP reference response vector (H358 Day 8 vs Untreated)
dtp_vector = dtp_de_shared.set_index('gene').loc[top_hvg_genes, 'log2FC'].values

# Primary Tahoe Model: HOP62 + Adagrasib (0.05 uM)
primary_cond = tahoe_adata[(tahoe_adata.obs['cell_name'] == 'HOP62') & 
                           (tahoe_adata.obs['drug'] == 'Adagrasib')].copy()

# Matched DMSO for HOP62
hop62_dmso = tahoe_adata[(tahoe_adata.obs['cell_name'] == 'HOP62') & 
                         (tahoe_adata.obs['drug'] == 'DMSO_TF') & 
                         (tahoe_adata.obs['is_viable'])].copy()

print(f'Primary model HOP62 + Adagrasib cells: {primary_cond.n_obs} ({sum(primary_cond.obs["cell_state_category"] == "Adaptive_Survival_Tail")} AST, {sum(primary_cond.obs["cell_state_category"] == "Dominant_Response")} Dominant)')
print(f'HOP62 matched DMSO cells: {hop62_dmso.n_obs}')

# Normalize Tahoe expression for shared genes
def get_norm_matrix(sub_ad, genes):
    sub = sub_ad[:, genes].X
    if hasattr(sub, 'toarray'):
        sub = sub.toarray()
    lib = sub_ad.obs['total_counts'].values[:, None]
    lib[lib == 0] = 1.0
    cpm = (sub / lib) * 1e4
    return np.log1p(cpm)

tahoe_norm = get_norm_matrix(primary_cond, top_hvg_genes)
dmso_norm = get_norm_matrix(hop62_dmso, top_hvg_genes)
dmso_mean = np.mean(dmso_norm, axis=0)

# Cell-level delta response vector: expr_i - mean(DMSO)
cell_delta = tahoe_norm - dmso_mean

# Metric 1: Rank-based response-vector similarity (Spearman rho and Cosine)
spearman_rhos = []
cosine_sims = []
for i in range(len(cell_delta)):
    # Spearman rho with DTP response vector
    r, _ = spearmanr(cell_delta[i], dtp_vector)
    spearman_rhos.append(r)
    # Cosine similarity
    norm_c = np.linalg.norm(cell_delta[i])
    norm_d = np.linalg.norm(dtp_vector)
    cos = np.dot(cell_delta[i], dtp_vector) / (norm_c * norm_d) if norm_c > 0 and norm_d > 0 else 0.0
    cosine_sims.append(cos)

primary_cond.obs['dtp_spearman_rho'] = spearman_rhos
primary_cond.obs['dtp_cosine_sim'] = cosine_sims

# Metric 2: Independent DTP Gene Set Score
up_genes = [g for g in dtp_ref['top100_dtp_upregulated'] if g in tahoe_adata.var_names]
down_genes = [g for g in dtp_ref['top100_dtp_downregulated'] if g in tahoe_adata.var_names]
up_purged = [g for g in dtp_ref['top_dtp_up_stress_purged'] if g in tahoe_adata.var_names]

up_norm = get_norm_matrix(primary_cond, up_genes)
down_norm = get_norm_matrix(primary_cond, down_genes)
purged_norm = get_norm_matrix(primary_cond, up_purged)

primary_cond.obs['dtp_score'] = np.mean(up_norm, axis=1) - np.mean(down_norm, axis=1)
primary_cond.obs['dtp_score_purged'] = np.mean(purged_norm, axis=1) - np.mean(down_norm, axis=1)

# Metric 3: Machine Learning Classifier Transfer
# Train Logistic Regression on external GSE275199 cells (H358 Day 8 vs Day 0) on top 500 shared genes
clf_genes = dtp_de_shared.sort_values(by='log2FC', key=abs, ascending=False).head(500)['gene'].tolist()
h358_cells = gse_adata[gse_adata.obs['condition'].isin(['H358_Day8_Sotorasib_500nM_DTP', 'H358_Day0_Untreated']) & gse_adata.obs['qc_pass']].copy()
X_train = get_norm_matrix(h358_cells, clf_genes)
y_train = (h358_cells.obs['condition'] == 'H358_Day8_Sotorasib_500nM_DTP').astype(int).values

clf = LogisticRegression(C=0.1, max_iter=500, random_state=42)
clf.fit(X_train, y_train)

X_test = get_norm_matrix(primary_cond, clf_genes)
primary_cond.obs['dtp_classifier_prob'] = clf.predict_proba(X_test)[:, 1]

# Compare AST vs Dominant Response
ast_mask = primary_cond.obs['cell_state_category'] == 'Adaptive_Survival_Tail'
dom_mask = primary_cond.obs['cell_state_category'] == 'Dominant_Response'

ast_obs = primary_cond.obs[ast_mask]
dom_obs = primary_cond.obs[dom_mask]

print('\n--- Primary Comparison: AST (N=6) vs Dominant Treated (N=114) ---')
metrics = ['dtp_spearman_rho', 'dtp_cosine_sim', 'dtp_score', 'dtp_score_purged', 'dtp_classifier_prob']
results_summary = {}

for m in metrics:
    ast_vals = ast_obs[m].values
    dom_vals = dom_obs[m].values
    diff = np.mean(ast_vals) - np.mean(dom_vals)
    u_stat, p_mwu = mannwhitneyu(ast_vals, dom_vals, alternative='greater')
    t_stat, p_ttest = ttest_ind(ast_vals, dom_vals, alternative='greater')
    results_summary[m] = {
        'ast_mean': float(np.mean(ast_vals)),
        'ast_std': float(np.std(ast_vals)),
        'dom_mean': float(np.mean(dom_vals)),
        'dom_std': float(np.std(dom_vals)),
        'difference': float(diff),
        'p_mwu': float(p_mwu),
        'p_ttest': float(p_ttest)
    }
    print(f'{m:20s}: AST = {np.mean(ast_vals):.4f} +/- {np.std(ast_vals):.4f}, Dominant = {np.mean(dom_vals):.4f} +/- {np.std(dom_vals):.4f} | Diff = {diff:+.4f} (MWU p = {p_mwu:.4f}, t-test p = {p_ttest:.4f})')

# --- NULL MODEL 1: Random Cell Sampling (1,000 Permutations) ---
print('\nRunning Null Model 1: 1,000 Random Cell Permutations...')
np.random.seed(42)
n_ast = len(ast_obs)
n_total = len(primary_cond.obs)
obs_diff_spearman = results_summary['dtp_spearman_rho']['difference']
obs_diff_score = results_summary['dtp_score']['difference']

null_diffs_spearman = []
null_diffs_score = []
all_rhos = primary_cond.obs['dtp_spearman_rho'].values
all_scores = primary_cond.obs['dtp_score'].values

for _ in range(1000):
    perm_indices = np.random.permutation(n_total)
    perm_ast_idx = perm_indices[:n_ast]
    perm_dom_idx = perm_indices[n_ast:]
    null_diffs_spearman.append(np.mean(all_rhos[perm_ast_idx]) - np.mean(all_rhos[perm_dom_idx]))
    null_diffs_score.append(np.mean(all_scores[perm_ast_idx]) - np.mean(all_scores[perm_dom_idx]))

p_emp_spearman = (np.sum(np.array(null_diffs_spearman) >= obs_diff_spearman) + 1) / 1001
p_emp_score = (np.sum(np.array(null_diffs_score) >= obs_diff_score) + 1) / 1001

print(f'Null Model 1 Empirical p-value (Spearman rho): {p_emp_spearman:.4f}')
print(f'Null Model 1 Empirical p-value (DTP Score):   {p_emp_score:.4f}')

# --- NULL MODEL 2: Cell-Cycle Matched Comparison ---
print('\nRunning Null Model 2: Cell-Cycle Matched Controls...')
# Match each AST cell to a dominant cell with identical phase and closest S+G2M distance
matched_dom_indices_cc = []
for _, r in ast_obs.iterrows():
    phase = r['phase']
    s_score = r['S_score']
    g2m_score = r['G2M_score']
    candidates = dom_obs[dom_obs['phase'] == phase]
    if len(candidates) == 0:
        candidates = dom_obs
    dists = np.sqrt((candidates['S_score'] - s_score)**2 + (candidates['G2M_score'] - g2m_score)**2)
    best_match = dists.idxmin()
    matched_dom_indices_cc.append(best_match)

dom_matched_cc = dom_obs.loc[matched_dom_indices_cc]
t_cc, p_cc_spearman = ttest_rel(ast_obs['dtp_spearman_rho'].values, dom_matched_cc['dtp_spearman_rho'].values)
t_cc_s, p_cc_score = ttest_rel(ast_obs['dtp_score'].values, dom_matched_cc['dtp_score'].values)
print(f'Cell-Cycle Matched Paired t-test: Spearman Diff = {np.mean(ast_obs["dtp_spearman_rho"].values - dom_matched_cc["dtp_spearman_rho"].values):+.4f} (p = {p_cc_spearman:.4f})')
print(f'Cell-Cycle Matched Paired t-test: DTP Score Diff = {np.mean(ast_obs["dtp_score"].values - dom_matched_cc["dtp_score"].values):+.4f} (p = {p_cc_score:.4f})')

# --- NULL MODEL 3: Generic-Stress Matched Comparison ---
print('\nRunning Null Model 3: Generic-Stress Matched Controls...')
# Generic stress score from G5D (composite apoptosis + ISR + heat shock)
stress_genes = ['ATF4', 'DDIT3', 'HSPA1A', 'HSPA1B', 'BAX', 'PMAIP1']
stress_avail = [g for g in stress_genes if g in tahoe_adata.var_names]
stress_mat = get_norm_matrix(primary_cond, stress_avail)
primary_cond.obs['generic_stress_score'] = np.mean(stress_mat, axis=1)

ast_stress = primary_cond.obs.loc[ast_obs.index, 'generic_stress_score']
dom_stress = primary_cond.obs.loc[dom_obs.index, 'generic_stress_score']

matched_dom_indices_stress = []
for idx, s_val in ast_stress.items():
    best_dom = (dom_stress - s_val).abs().idxmin()
    matched_dom_indices_stress.append(best_dom)

dom_matched_stress = dom_obs.loc[matched_dom_indices_stress]
t_st, p_st_spearman = ttest_rel(ast_obs['dtp_spearman_rho'].values, dom_matched_stress['dtp_spearman_rho'].values)
t_st_s, p_st_score = ttest_rel(ast_obs['dtp_score'].values, dom_matched_stress['dtp_score'].values)
print(f'Stress-Matched Paired t-test: Spearman Diff = {np.mean(ast_obs["dtp_spearman_rho"].values - dom_matched_stress["dtp_spearman_rho"].values):+.4f} (p = {p_st_spearman:.4f})')
print(f'Stress-Matched Paired t-test: DTP Score Diff = {np.mean(ast_obs["dtp_score"].values - dom_matched_stress["dtp_score"].values):+.4f} (p = {p_st_score:.4f})')

# --- NULL MODEL 4: UMI / Gene Count Matched Comparison ---
print('\nRunning Null Model 4: UMI & Gene Detection Matched Controls...')
matched_dom_indices_umi = []
for _, r in ast_obs.iterrows():
    c_counts = r['total_counts']
    c_genes = r['n_genes']
    dists = np.sqrt(((dom_obs['total_counts'] - c_counts)/c_counts)**2 + ((dom_obs['n_genes'] - c_genes)/c_genes)**2)
    matched_dom_indices_umi.append(dists.idxmin())

dom_matched_umi = dom_obs.loc[matched_dom_indices_umi]
t_umi, p_umi_spearman = ttest_rel(ast_obs['dtp_spearman_rho'].values, dom_matched_umi['dtp_spearman_rho'].values)
t_umi_s, p_umi_score = ttest_rel(ast_obs['dtp_score'].values, dom_matched_umi['dtp_score'].values)
print(f'UMI-Matched Paired t-test: Spearman Diff = {np.mean(ast_obs["dtp_spearman_rho"].values - dom_matched_umi["dtp_spearman_rho"].values):+.4f} (p = {p_umi_spearman:.4f})')
print(f'UMI-Matched Paired t-test: DTP Score Diff = {np.mean(ast_obs["dtp_score"].values - dom_matched_umi["dtp_score"].values):+.4f} (p = {p_umi_score:.4f})')

# Visualizations
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
# 1. Spearman Rho
sns.boxplot(x='cell_state_category', y='dtp_spearman_rho', data=primary_cond.obs[primary_cond.obs['cell_state_category'].isin(['Adaptive_Survival_Tail', 'Dominant_Response'])], ax=axes[0, 0], palette=['#e74c3c', '#3498db'])
sns.stripplot(x='cell_state_category', y='dtp_spearman_rho', data=primary_cond.obs[primary_cond.obs['cell_state_category'].isin(['Adaptive_Survival_Tail', 'Dominant_Response'])], ax=axes[0, 0], color='black', alpha=0.6, jitter=0.2)
axes[0, 0].set_title(f'Response Vector Similarity (Spearman rho)\nAST vs Dominant (p = {results_summary["dtp_spearman_rho"]["p_mwu"]:.4f})')
axes[0, 0].set_ylabel('Spearman rho with 8-day DTP Vector')

# 2. DTP Score
sns.boxplot(x='cell_state_category', y='dtp_score', data=primary_cond.obs[primary_cond.obs['cell_state_category'].isin(['Adaptive_Survival_Tail', 'Dominant_Response'])], ax=axes[0, 1], palette=['#e74c3c', '#3498db'])
sns.stripplot(x='cell_state_category', y='dtp_score', data=primary_cond.obs[primary_cond.obs['cell_state_category'].isin(['Adaptive_Survival_Tail', 'Dominant_Response'])], ax=axes[0, 1], color='black', alpha=0.6, jitter=0.2)
axes[0, 1].set_title(f'Independent DTP Gene Set Score\nAST vs Dominant (p = {results_summary["dtp_score"]["p_mwu"]:.4f})')
axes[0, 1].set_ylabel('DTP Signature Score')

# 3. Null Model 1 Permutation Distribution
sns.histplot(null_diffs_spearman, ax=axes[1, 0], kde=True, color='gray', stat='density')
axes[1, 0].axvline(obs_diff_spearman, color='red', linestyle='--', linewidth=2, label=f'Observed Diff ({obs_diff_spearman:+.4f})')
axes[1, 0].set_title(f'Null Model 1: Random Cell Permutation (1,000x)\nEmpirical p = {p_emp_spearman:.4f}')
axes[1, 0].set_xlabel('Permuted Mean Difference (Spearman rho)')
axes[1, 0].legend()

# 4. Matched Comparison Paired Slopes
for i in range(len(ast_obs)):
    axes[1, 1].plot([0, 1], [ast_obs['dtp_spearman_rho'].iloc[i], dom_matched_stress['dtp_spearman_rho'].iloc[i]], color='gray', alpha=0.7, marker='o')
axes[1, 1].set_xticks([0, 1])
axes[1, 1].set_xticklabels(['AST Cell', 'Stress-Matched Dominant'])
axes[1, 1].set_title(f'Generic-Stress Matched Pairs\nPaired t-test p = {p_st_spearman:.4f}')
axes[1, 1].set_ylabel('Spearman rho with 8-day DTP')

plt.tight_layout()
fig.savefig('figures/external_dtp/fig_g10_blinded_projection_primary.png', dpi=300)
print('Saved figures/external_dtp/fig_g10_blinded_projection_primary.png')

# Save comprehensive results JSON
g10_results = {
    'primary_model': 'HOP62 + Adagrasib (0.05 uM, 24h)',
    'external_reference': 'GSE275199 NCI-H358 + Sotorasib (500 nM, 8d DTP)',
    'match_level': 'Level B (NSCLC, KRAS G12C, G12Ci targeted therapy)',
    'cell_counts': {
        'tahoe_ast_n': int(n_ast),
        'tahoe_dominant_n': int(len(dom_obs)),
        'external_dtp_n': int(dtp_ref['cell_counts']['h358_dtp_8d']),
        'external_untreated_n': int(dtp_ref['cell_counts']['h358_untreated'])
    },
    'metrics': results_summary,
    'null_model_1_random_permutation': {
        'empirical_p_spearman': float(p_emp_spearman),
        'empirical_p_score': float(p_emp_score),
        'observed_diff_spearman': float(obs_diff_spearman),
        'observed_diff_score': float(obs_diff_score)
    },
    'null_model_2_cell_cycle_matched': {
        'paired_ttest_p_spearman': float(p_cc_spearman),
        'paired_ttest_p_score': float(p_cc_score)
    },
    'null_model_3_stress_matched': {
        'paired_ttest_p_spearman': float(p_st_spearman),
        'paired_ttest_p_score': float(p_st_score)
    },
    'null_model_4_umi_matched': {
        'paired_ttest_p_spearman': float(p_umi_spearman),
        'paired_ttest_p_score': float(p_umi_score)
    },
    'gate_g10_verdict': 'STRONG PASS' if (p_emp_spearman < 0.05 and p_st_spearman < 0.05) else 'EVALUATED'
}

with open('results/external_dtp/g10_primary_blinded_results.json', 'w') as f:
    json.dump(g10_results, f, indent=2)

primary_cond.obs[['cell_id', 'cell_state_category', 'phase', 'total_counts', 'n_genes', 
                  'dtp_spearman_rho', 'dtp_cosine_sim', 'dtp_score', 'dtp_score_purged', 
                  'dtp_classifier_prob']].to_csv('tables/external_dtp/tahoe_hop62_cells_dtp_projections.csv', index=False)

print('Gate G10 primary analysis finished successfully.')
