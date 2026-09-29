import scanpy as sc
import numpy as np
from scipy.stats import ttest_rel
import pandas as pd

tahoe_adata = sc.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
tahoe_adata.obs = tahoe_adata.obs.set_index('cell_id')
proj = pd.read_csv('tables/external_dtp/tahoe_hop62_cells_dtp_projections.csv').set_index('cell_id')

primary_cells = proj.index
primary_cond = tahoe_adata[primary_cells].copy()
primary_cond.obs = primary_cond.obs.join(proj[['cell_state_category', 'dtp_spearman_rho', 'dtp_score']])

ast_obs = primary_cond.obs[primary_cond.obs['cell_state_category'] == 'Adaptive_Survival_Tail']
dom_obs = primary_cond.obs[primary_cond.obs['cell_state_category'] == 'Dominant_Response']

stress_genes = ['ATF4', 'DDIT3', 'HSPA1A', 'HSPA1B', 'BAX', 'PMAIP1']
stress_avail = [g for g in stress_genes if g in tahoe_adata.var_names]

def get_norm_matrix(sub_ad, genes):
    sub = sub_ad[:, genes].X
    if hasattr(sub, 'toarray'):
        sub = sub.toarray()
    lib = sub_ad.obs['total_counts'].values[:, None]
    lib[lib == 0] = 1.0
    cpm = (sub / lib) * 1e4
    return np.log1p(cpm)

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
diffs = ast_obs['dtp_spearman_rho'].values - dom_matched_stress['dtp_spearman_rho'].values

print("Test Statistic (t):", round(float(t_st), 4))
print("Degrees of Freedom:", len(diffs) - 1)
print("Two-Sided p-value:", round(float(p_st_spearman), 6))
print("Mean Difference (AST - Matched Dom):", round(float(np.mean(diffs)), 6))
print("One-Sided p-value (AST > Dom):", round(float(1.0 - p_st_spearman/2.0), 6))
print("One-Sided p-value (AST < Dom):", round(float(p_st_spearman/2.0), 6))
for i in range(len(diffs)):
    a_val = ast_obs['dtp_spearman_rho'].values[i]
    d_val = dom_matched_stress['dtp_spearman_rho'].values[i]
    print("Pair " + str(i+1) + ": AST rho = " + str(round(a_val, 4)) + ", Dom rho = " + str(round(d_val, 4)) + ", Diff = " + str(round(diffs[i], 4)))
