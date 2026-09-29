import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import cdist
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs('code/benchmark', exist_ok=True)
os.makedirs('tables/benchmark', exist_ok=True)
os.makedirs('results/benchmark', exist_ok=True)
os.makedirs('figures/benchmark', exist_ok=True)

print("--- Running P1: Comparator Metric Benchmark (ARTI vs C0–C5) ---")

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

# 5 PCs from PCA for geometric distances
from sklearn.decomposition import PCA
pca = PCA(n_components=5, random_state=42)
pca_X = pca.fit_transform(norm_X)

astr_df = pd.read_csv('tables/astr/astr_all_conditions.csv')

# Precompute DMSO data per cell line and plate
dmso_profiles = {}
for (cname, plate), grp in qc_pass_obs[qc_pass_obs['drug'] == 'DMSO_TF'].groupby(['cell_name', 'plate'], observed=True):
    idx = grp['cell_idx'].values
    dmso_profiles[(cname, plate)] = {
        'norm_mean': norm_X[idx].mean(axis=0),
        'pca_X': pca_X[idx],
        'pca_mean': pca_X[idx].mean(axis=0),
        'count': len(idx)
    }

records = []

for idx, r in astr_df.iterrows():
    cname = r['cell_name']
    plate = r['plate']
    cond_id = r['condition_id']
    arti = r['astr_score']
    drug = r['drug']
    
    dmso = dmso_profiles.get((cname, plate))
    if dmso is None or dmso['count'] < 5:
        continue
        
    cond_cells = qc_pass_obs[(qc_pass_obs['cell_name'] == cname) & (qc_pass_obs['drug'] == drug) & (qc_pass_obs['plate'] == plate)]
    if len(cond_cells) < 5:
        continue
        
    c_idx = cond_cells['cell_idx'].values
    X_cond = norm_X[c_idx]
    P_cond = pca_X[c_idx]
    P_dmso = dmso['pca_X']
    
    # C0: Pseudobulk L2 distance (AR1)
    cond_mean = X_cond.mean(axis=0)
    c0_pseudobulk_l2 = float(np.linalg.norm(cond_mean - dmso['norm_mean']))
    
    # C1: Energy Distance (E-distance) in PCA space
    p_c_sub = P_cond if len(P_cond) <= 200 else P_cond[np.random.choice(len(P_cond), 200, replace=False)]
    p_d_sub = P_dmso if len(P_dmso) <= 200 else P_dmso[np.random.choice(len(P_dmso), 200, replace=False)]
    
    dist_xy = cdist(p_c_sub, p_d_sub, metric='euclidean')
    dist_xx = cdist(p_c_sub, p_c_sub, metric='euclidean')
    dist_yy = cdist(p_d_sub, p_d_sub, metric='euclidean')
    
    e_distance = float(2.0 * dist_xy.mean() - dist_xx.mean() - dist_yy.mean())
    
    # C2: Perturbation-Response Score (PS)
    delta_vec = P_cond.mean(axis=0) - dmso['pca_mean']
    norm_delta = np.linalg.norm(delta_vec)
    if norm_delta > 0:
        unit_vec = delta_vec / norm_delta
        projections = np.dot(P_cond - dmso['pca_mean'], unit_vec)
        dmso_projections = np.dot(P_dmso - dmso['pca_mean'], unit_vec)
        penetrance_threshold = np.percentile(dmso_projections, 95)
        response_penetrance_fraction = float(np.mean(projections > penetrance_threshold))
        mean_projection_magnitude = float(np.mean(projections))
        ps_fraction = response_penetrance_fraction
        ps_mean_magnitude = mean_projection_magnitude
    else:
        response_penetrance_fraction = 0.0
        mean_projection_magnitude = 0.0
        ps_fraction = 0.0
        ps_mean_magnitude = 0.0
        
    # C3: Simple within-condition dispersion
    cond_centroid = P_cond.mean(axis=0)
    c3_simple_dispersion = float(np.mean(np.linalg.norm(P_cond - cond_centroid, axis=1)))
    
    # C4: Transcriptomic Entropy
    cell_sums = X_cond.sum(axis=1, keepdims=True)
    cell_sums[cell_sums == 0] = 1.0
    p_mat = X_cond / cell_sums
    entropy_per_cell = -np.sum(np.where(p_mat > 0, p_mat * np.log(p_mat + 1e-12), 0), axis=1)
    c4_mean_entropy = float(np.mean(entropy_per_cell))
    
    # C5: Mean distance to control centroid
    c5_dist_to_ctrl_centroid = float(np.mean(np.linalg.norm(P_cond - dmso['pca_mean'], axis=1)))
    
    records.append({
        'condition_id': cond_id,
        'cell_name': cname,
        'drug': r['drug'],
        'plate': plate,
        'n_qc_pass': len(cond_cells),
        'arti': arti,
        'c0_pseudobulk_l2': c0_pseudobulk_l2,
        'c1_e_distance': e_distance,
        'c2_response_penetrance_fraction': response_penetrance_fraction,
        'c2_mean_projection_magnitude': mean_projection_magnitude,
        'c2_ps_fraction': ps_fraction,
        'c2_ps_mean_magnitude': ps_mean_magnitude,
        'c3_simple_dispersion': c3_simple_dispersion,
        'c4_transcriptomic_entropy': c4_mean_entropy,
        'c5_dist_to_control_centroid': c5_dist_to_ctrl_centroid
    })

bench_df = pd.DataFrame(records)

comparators = [
    'c0_pseudobulk_l2', 'c1_e_distance', 'c2_response_penetrance_fraction', 
    'c2_mean_projection_magnitude', 'c3_simple_dispersion', 
    'c4_transcriptomic_entropy', 'c5_dist_to_control_centroid'
]

bench_df = bench_df.dropna(subset=['arti'] + comparators).reset_index(drop=True)
bench_df.to_csv('tables/benchmark/comparator_metric_benchmark.csv', index=False)
print("Benchmarked", len(bench_df), "valid conditions.")

summary_stats = {}
for comp in comparators:
    pr, p_val = pearsonr(bench_df['arti'], bench_df[comp])
    sr, s_val = spearmanr(bench_df['arti'], bench_df[comp])
    r2 = pr ** 2
    summary_stats[comp] = {
        'pearson_r': float(pr),
        'pearson_p': float(p_val),
        'spearman_rho': float(sr),
        'spearman_p': float(s_val),
        'r_squared': float(r2)
    }
    print(f"ARTI vs {comp:32s}: Pearson r = {pr:+.4f} (p = {p_val:.2e}), Spearman rho = {sr:+.4f}, R2 = {r2:.4f}")

# Backwards compatibility aliases
summary_stats['c2_ps_fraction'] = summary_stats['c2_response_penetrance_fraction']
summary_stats['c2_ps_mean_magnitude'] = summary_stats['c2_mean_projection_magnitude']

# Multivariate OLS model
X_reg = bench_df[comparators]
X_reg = sm.add_constant(X_reg)
y_reg = bench_df['arti']
model = sm.OLS(y_reg, X_reg).fit()

residual_variance_fraction = float(1.0 - model.rsquared)
print(f"\nMultivariate Regression R2 (ARTI ~ all comparators): {model.rsquared:.4f}")
print(f"Unexplained / Independent Residual Variance: {residual_variance_fraction * 100:.2f}%")

bench_df['arti_residual'] = model.resid
bench_df.to_csv('tables/benchmark/comparator_metric_benchmark.csv', index=False)
bench_df.to_csv('tables/benchmark/final_metric_comparison.csv', index=False)

# E-distance empirical decomposition (P13.7)
r_e_c0, p_e_c0 = pearsonr(bench_df['c1_e_distance'], bench_df['c0_pseudobulk_l2'])
r_e_c5, p_e_c5 = pearsonr(bench_df['c1_e_distance'], bench_df['c5_dist_to_control_centroid'])
r_e_c3, p_e_c3 = pearsonr(bench_df['c1_e_distance'], bench_df['c3_simple_dispersion'])
r_e_arti, p_e_arti = pearsonr(bench_df['c1_e_distance'], bench_df['arti'])

print(f"\n--- P13.7 E-distance Empirical Decomposition ---")
print(f"  E-distance vs C0 (Pseudobulk L2): r = {r_e_c0:+.4f} (p = {p_e_c0:.2e})")
print(f"  E-distance vs C5 (Ctrl Centroid): r = {r_e_c5:+.4f} (p = {p_e_c5:.2e})")
print(f"  E-distance vs C3 (Within Disp):   r = {r_e_c3:+.4f} (p = {p_e_c3:.2e})")
print(f"  E-distance vs ARTI:               r = {r_e_arti:+.4f} (p = {p_e_arti:.2e})")

# Nonlinear Redundancy Test via 5-fold Cross-Validation (P13.6)
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import KFold, cross_val_score

X_mat = bench_df[comparators].values
y_vec = bench_df['arti'].values
kf = KFold(n_splits=5, shuffle=True, random_state=42)
rf_model = RandomForestRegressor(n_estimators=100, max_depth=4, random_state=42)
gb_model = HistGradientBoostingRegressor(max_depth=3, random_state=42)

cv_r2_rf = float(np.mean(cross_val_score(rf_model, X_mat, y_vec, cv=kf, scoring='r2')))
cv_r2_gb = float(np.mean(cross_val_score(gb_model, X_mat, y_vec, cv=kf, scoring='r2')))

print(f"\n--- P13.6 Nonlinear Redundancy CV R2 ---")
print(f"  Random Forest 5-fold CV R2:      {cv_r2_rf:.4f}")
print(f"  Gradient Boosting 5-fold CV R2:  {cv_r2_gb:.4f}")


results_summary = {
    'conditions_evaluated': len(bench_df),
    'bivariate_correlations': summary_stats,
    'multivariate_model': {
        'r_squared': float(model.rsquared),
        'adj_r_squared': float(model.rsquared_adj),
        'f_statistic': float(model.fvalue),
        'f_pvalue': float(model.f_pvalue),
        'residual_orthogonal_variance_fraction': residual_variance_fraction
    },
    'e_distance_decomposition': {
        'r_vs_c0_l2': float(r_e_c0),
        'r_vs_c5_ctrl_centroid': float(r_e_c5),
        'r_vs_c3_dispersion': float(r_e_c3),
        'r_vs_arti': float(r_e_arti)
    },
    'nonlinear_redundancy_cv': {
        'random_forest_cv_r2': cv_r2_rf,
        'gradient_boosting_cv_r2': cv_r2_gb
    }
}

with open('results/benchmark/comparator_benchmark_summary.json', 'w') as f:
    json.dump(results_summary, f, indent=2)

# Figure: Multi-panel scatter plot
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('Project A5: ARTI vs Existing Single-Cell Perturbation Comparators (N = 263)', fontsize=15, y=0.98)

plot_items = [
    ('c0_pseudobulk_l2', 'C0: Pseudobulk L2 Distance', axes[0, 0]),
    ('c1_e_distance', 'C1: Energy Distance (scPerturb)', axes[0, 1]),
    ('c2_response_penetrance_fraction', 'C2: Response Penetrance Fraction', axes[0, 2]),
    ('c3_simple_dispersion', 'C3: Within-Condition Dispersion', axes[1, 0]),
    ('c4_transcriptomic_entropy', 'C4: Transcriptomic Entropy', axes[1, 1]),
    ('c5_dist_to_control_centroid', 'C5: Dist to Control Centroid', axes[1, 2])
]

for col, title, ax in plot_items:
    r_val = summary_stats[col]['pearson_r']
    r2_val = summary_stats[col]['r_squared']
    sns.regplot(data=bench_df, x=col, y='arti', ax=ax,
                scatter_kws={'alpha': 0.6, 'color': '#2b5c8f', 's': 30},
                line_kws={'color': '#d95f02', 'linewidth': 2})
    ax.set_title(f"{title}\n(Pearson r = {r_val:+.3f}, R² = {r2_val:.3f})", fontsize=11)
    ax.set_ylabel('ARTI (Acute Response Tail Index)', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('figures/benchmark/fig_comparator_correlations.png', dpi=300)
plt.close()
print("Saved figure to figures/benchmark/fig_comparator_correlations.png")
