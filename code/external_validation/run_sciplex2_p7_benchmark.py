import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, ttest_rel, ttest_ind
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import statsmodels.api as sm
import anndata as ad

print("=== Running P7 External Acute Perturbation Benchmark on SciPlex2 ===")
t0 = time.time()

# 1. Load SciPlex2 AnnData
adata_path = 'data_processed/sciplex2_adata.h5ad'
adata = ad.read_h5ad(adata_path)
print(f"Loaded {adata.shape[0]} cells, {adata.shape[1]} genes.")

obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

# QC filtering: non-null drug & dose, UMI >= 500
valid_mask = obs['top_oligo'].notna() & (obs['n.umi'] >= 500)
qc_obs = obs[valid_mask].copy()
print(f"QC pass cells: {len(qc_obs)} / {len(obs)} ({len(qc_obs)/len(obs)*100:.1f}%)")

# Standardize well 'E09' to 'AE09' in Plate A
qc_obs['well'] = qc_obs['well'].astype(str).replace({'E09': 'AE09'})
qc_obs['plate'] = qc_obs['well'].str[0]

# 2. Top 1000 HVGs (frozen definition v1)
print("Selecting top 1000 HVGs and computing 5-PC PCA...")
sub_qc_idx = qc_obs['cell_idx'].values
X_qc = adata.X[sub_qc_idx]

gene_means = np.array(X_qc.mean(axis=0)).flatten()
gene_sq_means = np.array(X_qc.power(2).mean(axis=0)).flatten()
gene_vars = gene_sq_means - np.power(gene_means, 2)

top_1000_genes = np.argsort(gene_vars)[-1000:]
X_hvg = adata.X[:, top_1000_genes].toarray()

lib_sizes = X_hvg.sum(axis=1, keepdims=True)
lib_sizes[lib_sizes == 0] = 1.0
norm_X = np.log1p((X_hvg / lib_sizes) * 10000.0)

pca = PCA(n_components=5, random_state=42)
pca.fit(norm_X[sub_qc_idx])
pca_X = pca.transform(norm_X)
print(f"PCA explained variance ratio (5 PCs): {pca.explained_variance_ratio_.sum()*100:.2f}%")

# 3. Precompute vehicle control baselines per drug and plate
drugs = ['Nutlin', 'Dex', 'SAHA', 'BMS']
vehicle_profiles = {}

for drug in drugs:
    v_cells = qc_obs[(qc_obs['drug'] == drug) & (qc_obs['dose'] == '0')]
    v_idx = v_cells['cell_idx'].values
    v_pca = pca_X[v_idx]
    v_centroid = np.median(v_pca, axis=0, keepdims=True)
    v_dists = np.linalg.norm(v_pca - v_centroid, axis=1)
    v_disp = float(np.mean(v_dists))
    
    v_A_cells = v_cells[v_cells['plate'] == 'A']
    v_A_idx = v_A_cells['cell_idx'].values
    v_A_pca = pca_X[v_A_idx]
    v_A_disp = float(np.mean(np.linalg.norm(v_A_pca - np.median(v_A_pca, axis=0, keepdims=True), axis=1))) if len(v_A_idx) >= 10 else v_disp
    
    v_B_cells = v_cells[v_cells['plate'] == 'B']
    v_B_idx = v_B_cells['cell_idx'].values
    v_B_pca = pca_X[v_B_idx]
    v_B_disp = float(np.mean(np.linalg.norm(v_B_pca - np.median(v_B_pca, axis=0, keepdims=True), axis=1))) if len(v_B_idx) >= 10 else v_disp
    
    vehicle_profiles[drug] = {
        'all_cells': v_cells,
        'all_idx': v_idx,
        'all_pca': v_pca,
        'all_norm': norm_X[v_idx],
        'all_disp': v_disp,
        'plate_A_disp': v_A_disp,
        'plate_B_disp': v_B_disp,
        'norm_mean': norm_X[v_idx].mean(axis=0),
        'pca_mean': v_pca.mean(axis=0),
        'pca_centroid': v_centroid[0]
    }
    print(f"Vehicle baseline {drug:>6}: N={len(v_cells)}, overall disp={v_disp:.4f}, Plate A={v_A_disp:.4f}, Plate B={v_B_disp:.4f}")

# 4. Compute condition metrics for all wells
records = []
for (drug, dose, well), grp in qc_obs.groupby(['drug', 'dose', 'well']):
    n_cells = len(grp)
    if n_cells < 15:
        continue
    c_idx = grp['cell_idx'].values
    plate = grp['plate'].iloc[0]
    dose_num = grp['dose_numeric'].iloc[0]
    
    v_prof = vehicle_profiles[drug]
    matched_v_disp = v_prof['plate_A_disp'] if plate == 'A' else v_prof['plate_B_disp']
    
    P_cond = pca_X[c_idx]
    cond_median = np.median(P_cond, axis=0, keepdims=True)
    dists = np.linalg.norm(P_cond - cond_median, axis=1)
    cond_disp = float(np.mean(dists))
    
    delta_disp = cond_disp - matched_v_disp
    arti = delta_disp / 2.0
    
    v_dists = np.linalg.norm(v_prof['all_pca'] - v_prof['pca_centroid'], axis=1)
    v_th_95 = np.percentile(v_dists, 95.0)
    tail_fraction = float(np.mean(dists > v_th_95))
    
    X_cond = norm_X[c_idx]
    cond_mean = X_cond.mean(axis=0)
    c0_pseudobulk_l2 = float(np.linalg.norm(cond_mean - v_prof['norm_mean']))
    
    P_v = v_prof['all_pca']
    sub_c = P_cond if len(P_cond) <= 150 else P_cond[np.random.choice(len(P_cond), 150, replace=False)]
    sub_v = P_v if len(P_v) <= 150 else P_v[np.random.choice(len(P_v), 150, replace=False)]
    d_xy = cdist(sub_c, sub_v, metric='euclidean')
    d_xx = cdist(sub_c, sub_c, metric='euclidean')
    d_yy = cdist(sub_v, sub_v, metric='euclidean')
    c1_energy_dist = float(2.0 * d_xy.mean() - d_xx.mean() - d_yy.mean())
    
    delta_vec = P_cond.mean(axis=0) - v_prof['pca_mean']
    norm_delta = np.linalg.norm(delta_vec)
    if norm_delta > 0:
        u_vec = delta_vec / norm_delta
        proj_c = np.dot(P_cond - v_prof['pca_mean'], u_vec)
        proj_v = np.dot(P_v - v_prof['pca_mean'], u_vec)
        ps_th = np.percentile(proj_v, 95.0)
        c2_ps_fraction = float(np.mean(proj_c > ps_th))
        c2_ps_magnitude = float(np.mean(proj_c))
    else:
        c2_ps_fraction = 0.0
        c2_ps_magnitude = 0.0
        
    c3_simple_disp = cond_disp
    
    records.append({
        'condition_id': f"{drug}_{dose}_{well}",
        'drug': drug,
        'dose': str(dose),
        'dose_numeric': dose_num,
        'well': well,
        'plate': plate,
        'n_cells': n_cells,
        'dispersion': cond_disp,
        'matched_v_dispersion': matched_v_disp,
        'delta_dispersion': delta_disp,
        'arti': arti,
        'tail_fraction': tail_fraction,
        'c0_pseudobulk_l2': c0_pseudobulk_l2,
        'c1_energy_dist': c1_energy_dist,
        'c2_ps_fraction': c2_ps_fraction,
        'c2_ps_magnitude': c2_ps_magnitude,
        'c3_simple_dispersion': c3_simple_disp
    })

df_metrics = pd.DataFrame(records)
os.makedirs('tables/benchmark', exist_ok=True)
df_metrics.to_csv('tables/benchmark/sciplex2_condition_metrics.csv', index=False)
print(f"Computed metrics for {len(df_metrics)} well conditions. Saved tables/benchmark/sciplex2_condition_metrics.csv.")

# Filter to active perturbations (dose > 0)
df_active = df_metrics[df_metrics['dose_numeric'] > 0].copy()
print(f"Active perturbation conditions: {len(df_active)}")

# 5. Replicate Reproducibility Analysis (Cross-Plate & Cross-Well)
# For each drug and dose, compare Plate A vs Plate B mean ARTI
rep_pairs = []
for (drug, dose), grp in df_active.groupby(['drug', 'dose']):
    grp_A = grp[grp['plate'] == 'A']
    grp_B = grp[grp['plate'] == 'B']
    if len(grp_A) > 0 and len(grp_B) > 0:
        rep_pairs.append({
            'drug': drug,
            'dose': dose,
            'dose_numeric': grp['dose_numeric'].iloc[0],
            'arti_plate_A': grp_A['arti'].mean(),
            'arti_plate_B': grp_B['arti'].mean(),
            'c0_plate_A': grp_A['c0_pseudobulk_l2'].mean(),
            'c0_plate_B': grp_B['c0_pseudobulk_l2'].mean(),
            'n_wells_A': len(grp_A),
            'n_wells_B': len(grp_B)
        })

df_reps = pd.DataFrame(rep_pairs)
df_reps.to_csv('tables/benchmark/sciplex2_replicate_pairs.csv', index=False)

r_p, p_p = pearsonr(df_reps['arti_plate_A'], df_reps['arti_plate_B'])
r_s, p_s = spearmanr(df_reps['arti_plate_A'], df_reps['arti_plate_B'])

# ICC(2,1) and ICC(3,1)
Y = df_reps[['arti_plate_A', 'arti_plate_B']].values
n_items = len(Y)
mean_per_item = Y.mean(axis=1)
grand_mean = Y.mean()
ss_total = np.sum((Y - grand_mean)**2)
ss_between_items = 2.0 * np.sum((mean_per_item - grand_mean)**2)
ss_within_items = ss_total - ss_between_items
mean_per_judge = Y.mean(axis=0)
ss_between_judges = n_items * np.sum((mean_per_judge - grand_mean)**2)
ss_error = ss_within_items - ss_between_judges
ms_items = ss_between_items / (n_items - 1)
ms_judges = ss_between_judges / 1
ms_error = ss_error / (n_items - 1)
icc_21 = float((ms_items - ms_error) / (ms_items + ms_error + 2.0 * (ms_judges - ms_error) / n_items))
icc_31 = float((ms_items - ms_error) / (ms_items + ms_error))

print("\n--- SciPlex2 Replicate Reproducibility Results ---")
print(f"Drug-dose pairs evaluated (Plate A vs B): {len(df_reps)}")
print(f"Pearson r = {r_p:+.4f} (p = {p_p:.4e})")
print(f"Spearman rho = {r_s:+.4f} (p = {p_s:.4e})")
print(f"ICC(2,1) = {icc_21:+.4f}, ICC(3,1) = {icc_31:+.4f}")

# 6. Dose-Response Architecture Analysis
dose_summary = []
for drug, grp in df_metrics.groupby('drug'):
    # Mean ARTI and C0 per dose
    dose_means = grp.groupby('dose_numeric').agg({
        'arti': ['mean', 'std', 'count'],
        'c0_pseudobulk_l2': 'mean',
        'c1_energy_dist': 'mean',
        'c2_ps_fraction': 'mean'
    }).reset_index()
    dose_means.columns = ['dose_numeric', 'arti_mean', 'arti_std', 'n_wells', 'c0_mean', 'c1_mean', 'c2_mean']
    dose_means['drug'] = drug
    dose_summary.append(dose_means)
    
    # Dose correlation within drug
    active_grp = grp[grp['dose_numeric'] > 0]
    r_dose, p_dose = spearmanr(active_grp['dose_numeric'], active_grp['arti'])
    print(f"Dose response {drug:>6}: Spearman rho(dose, ARTI) = {r_dose:+.4f} (p = {p_dose:.4e})")

df_dose_arch = pd.concat(dose_summary, ignore_index=True)
df_dose_arch.to_csv('tables/benchmark/sciplex2_dose_architecture.csv', index=False)

# Overall dose correlation across all active conditions
r_dose_all, p_dose_all = spearmanr(df_active['dose_numeric'], df_active['arti'])
print(f"Overall dose-ARTI correlation: Spearman rho = {r_dose_all:+.4f} (p = {p_dose_all:.4e})")

# 7. Comparator Decoupling Analysis
corr_c0, p_c0 = pearsonr(df_active['arti'], df_active['c0_pseudobulk_l2'])
corr_c1, p_c1 = pearsonr(df_active['arti'], df_active['c1_energy_dist'])
corr_c2_frac, p_c2_frac = pearsonr(df_active['arti'], df_active['c2_ps_fraction'])
corr_c2_mag, p_c2_mag = pearsonr(df_active['arti'], df_active['c2_ps_magnitude'])
corr_c3, p_c3 = pearsonr(df_active['arti'], df_active['c3_simple_dispersion'])

# OLS Regression: ARTI ~ C0 + C1 + C2_mag + C3
X_ols = df_active[['c0_pseudobulk_l2', 'c1_energy_dist', 'c2_ps_magnitude', 'c3_simple_dispersion']]
X_ols_const = sm.add_constant(X_ols)
y_ols = df_active['arti']
model = sm.OLS(y_ols, X_ols_const).fit()
r2_multi = model.rsquared
residual_var = (1.0 - r2_multi) * 100.0

print("\n--- SciPlex2 Comparator Decoupling Results ---")
print(f"r(ARTI, C0 Pseudobulk L2) = {corr_c0:+.4f} (R2 = {corr_c0**2:.4f}, p = {p_c0:.4e})")
print(f"r(ARTI, C1 Energy Dist)   = {corr_c1:+.4f} (R2 = {corr_c1**2:.4f}, p = {p_c1:.4e})")
print(f"r(ARTI, C2 PS Magnitude)  = {corr_c2_mag:+.4f} (R2 = {corr_c2_mag**2:.4f}, p = {p_c2_mag:.4e})")
print(f"r(ARTI, C2 PS Fraction)   = {corr_c2_frac:+.4f} (R2 = {corr_c2_frac**2:.4f}, p = {p_c2_frac:.4e})")
print(f"r(ARTI, C3 Simple Disp)   = {corr_c3:+.4f} (R2 = {corr_c3**2:.4f}, p = {p_c3:.4e})")
print(f"Multivariate OLS R2 on all comparators: {r2_multi:.4f}")
print(f"ARTI Independent Residual Variance: {residual_var:.2f}%")

# 8. Vehicle Null False-Tail Rate Calibration
df_vehicle = df_metrics[df_metrics['dose_numeric'] == 0].copy()
v_artis = df_vehicle['arti'].values
false_tail_rate = float(np.mean(df_vehicle['arti'] > 0.05)) # threshold of significant tail expansion
mean_v_arti = float(np.mean(v_artis))
std_v_arti = float(np.std(v_artis))

print("\n--- SciPlex2 Vehicle Null Calibration ---")
print(f"Vehicle wells evaluated: {len(df_vehicle)}")
print(f"Vehicle mean ARTI: {mean_v_arti:+.4f} +/- {std_v_arti:.4f}")
print(f"Vehicle false positive tail rate (ARTI > 0.05): {false_tail_rate*100:.2f}% ({np.sum(df_vehicle['arti'] > 0.05)} / {len(df_vehicle)})")

# Save comparator summary table
comp_summary = pd.DataFrame([
    {'comparator': 'C0_Pseudobulk_L2', 'pearson_r': corr_c0, 'r_squared': corr_c0**2, 'p_value': p_c0},
    {'comparator': 'C1_Energy_Distance', 'pearson_r': corr_c1, 'r_squared': corr_c1**2, 'p_value': p_c1},
    {'comparator': 'C2_PS_Magnitude', 'pearson_r': corr_c2_mag, 'r_squared': corr_c2_mag**2, 'p_value': p_c2_mag},
    {'comparator': 'C2_PS_Fraction', 'pearson_r': corr_c2_frac, 'r_squared': corr_c2_frac**2, 'p_value': p_c2_frac},
    {'comparator': 'C3_Simple_Dispersion', 'pearson_r': corr_c3, 'r_squared': corr_c3**2, 'p_value': p_c3},
    {'comparator': 'Multivariate_OLS_All', 'pearson_r': np.sqrt(r2_multi), 'r_squared': r2_multi, 'p_value': model.f_pvalue}
])
comp_summary.to_csv('tables/benchmark/sciplex2_comparator_benchmark.csv', index=False)

t1 = time.time()
print(f"\nCompleted SciPlex2 P7 Benchmark in {t1-t0:.1f}s.")
