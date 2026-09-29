import os
import sys
import time
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
import anndata as ad

print("=== Running P7 External Acute Perturbation Benchmark on Aissa et al. 2021 ===")
t0 = time.time()

adata_path = 'data_raw/external_acute/AissaBenevolenskaya2021.h5ad'
adata = ad.read_h5ad(adata_path)
print(f"Loaded {adata.shape[0]} cells, {adata.shape[1]} genes.")

obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

# QC filtering: percent_mito < 15, ncounts >= 800, ngenes >= 300
qc_mask = (obs['percent_mito'] < 15.0) & (obs['ncounts'] >= 800) & (obs['ngenes'] >= 300)
qc_obs = obs[qc_mask].copy()
print(f"QC pass cells: {len(qc_obs)} / {len(obs)} ({len(qc_obs)/len(obs)*100:.1f}%)")
print(qc_obs['perturbation'].value_counts())

sub_idx = qc_obs['cell_idx'].values
np.random.seed(42)
sample_n = min(5000, len(sub_idx))
sample_idx = np.random.choice(sub_idx, sample_n, replace=False)
X_sub = adata.X[sample_idx]

gene_means = np.array(X_sub.mean(axis=0)).flatten()
gene_sq = np.array(X_sub.power(2).mean(axis=0)).flatten()
gene_vars = gene_sq - np.power(gene_means, 2)
top_1000_genes = np.argsort(gene_vars)[-1000:]

# Extract top 1000 genes for all QC cells
X_hvg = adata.X[sub_idx][:, top_1000_genes].toarray()
lib_sizes = X_hvg.sum(axis=1, keepdims=True)
lib_sizes[lib_sizes == 0] = 1.0
norm_X = np.log1p((X_hvg / lib_sizes) * 10000.0)

# Map sub_idx to row in norm_X
cell_to_norm_idx = {c: i for i, c in enumerate(sub_idx)}

# Fit 5-PC PCA on QC cells
pca = PCA(n_components=5, random_state=42)
pca_X = pca.fit_transform(norm_X)
print(f"PCA explained variance ratio (5 PCs): {pca.explained_variance_ratio_.sum()*100:.2f}%")

# Partition each condition into 5 independent non-overlapping replicates of N=250 cells
n_reps = 5
rep_size = 250
perturbations = ['control', 'crizotinib', 'osimertinib', 'osimertinib+crizotinib']

# Control baselines
ctrl_cells = qc_obs[qc_obs['perturbation'] == 'control'].copy()
ctrl_indices = [cell_to_norm_idx[c] for c in ctrl_cells['cell_idx'].values]
np.random.shuffle(ctrl_indices)

ctrl_reps = []
for r in range(n_reps):
    r_idx = ctrl_indices[r*rep_size : (r+1)*rep_size]
    P_ctrl = pca_X[r_idx]
    c_med = np.median(P_ctrl, axis=0, keepdims=True)
    disp = float(np.mean(np.linalg.norm(P_ctrl - c_med, axis=1)))
    ctrl_reps.append({
        'rep_id': r,
        'norm_idx': r_idx,
        'dispersion': disp,
        'pca_mean': P_ctrl.mean(axis=0),
        'norm_mean': norm_X[r_idx].mean(axis=0)
    })
mean_ctrl_disp = float(np.mean([cr['dispersion'] for cr in ctrl_reps]))

records = []
for pert in ['crizotinib', 'osimertinib', 'osimertinib+crizotinib']:
    p_cells = qc_obs[qc_obs['perturbation'] == pert].copy()
    p_indices = [cell_to_norm_idx[c] for c in p_cells['cell_idx'].values]
    np.random.shuffle(p_indices)
    
    for r in range(n_reps):
        r_idx = p_indices[r*rep_size : (r+1)*rep_size]
        P_cond = pca_X[r_idx]
        cond_median = np.median(P_cond, axis=0, keepdims=True)
        dists = np.linalg.norm(P_cond - cond_median, axis=1)
        disp = float(np.mean(dists))
        
        # Matched control replicate
        c_prof = ctrl_reps[r]
        delta_disp = disp - c_prof['dispersion']
        arti = delta_disp / 2.0
        
        # C0: Pseudobulk L2
        X_cond = norm_X[r_idx]
        c0 = float(np.linalg.norm(X_cond.mean(axis=0) - c_prof['norm_mean']))
        
        # C1: Energy distance
        P_ctrl = pca_X[c_prof['norm_idx']]
        sub_c = P_cond[:150]
        sub_v = P_ctrl[:150]
        d_xy = cdist(sub_c, sub_v, metric='euclidean')
        d_xx = cdist(sub_c, sub_c, metric='euclidean')
        d_yy = cdist(sub_v, sub_v, metric='euclidean')
        c1 = float(2.0 * d_xy.mean() - d_xx.mean() - d_yy.mean())
        
        # C2: PS
        delta_vec = P_cond.mean(axis=0) - c_prof['pca_mean']
        norm_delta = np.linalg.norm(delta_vec)
        if norm_delta > 0:
            u_vec = delta_vec / norm_delta
            proj_c = np.dot(P_cond - c_prof['pca_mean'], u_vec)
            proj_v = np.dot(P_ctrl - c_prof['pca_mean'], u_vec)
            ps_th = np.percentile(proj_v, 95.0)
            c2_frac = float(np.mean(proj_c > ps_th))
            c2_mag = float(np.mean(proj_c))
        else:
            c2_frac = 0.0
            c2_mag = 0.0
            
        records.append({
            'perturbation': pert,
            'replicate_id': r,
            'n_cells': rep_size,
            'dispersion': disp,
            'ctrl_dispersion': c_prof['dispersion'],
            'delta_dispersion': delta_disp,
            'arti': arti,
            'c0_pseudobulk_l2': c0,
            'c1_energy_dist': c1,
            'c2_ps_magnitude': c2_mag,
            'c2_ps_fraction': c2_frac
        })

df_aissa = pd.DataFrame(records)
df_aissa.to_csv('tables/benchmark/aissa_external_benchmark.csv', index=False)
print("Saved tables/benchmark/aissa_external_benchmark.csv.")

# Summary statistics per perturbation
print("\n--- Aissa et al. 2021 External Benchmark Summary ---")
for pert, grp in df_aissa.groupby('perturbation'):
    print(f"Perturbation {pert:>25}: ARTI = {grp['arti'].mean():+.4f} +/- {grp['arti'].std():.4f}, C0 = {grp['c0_pseudobulk_l2'].mean():.4f}, C1 = {grp['c1_energy_dist'].mean():.4f}")

# Cross-replicate correlation across perturbations (rep 0,1 vs rep 2,3)
split_1 = df_aissa[df_aissa['replicate_id'].isin([0, 1])].groupby('perturbation')['arti'].mean()
split_2 = df_aissa[df_aissa['replicate_id'].isin([2, 3])].groupby('perturbation')['arti'].mean()
r_rep, p_rep = pearsonr(split_1, split_2)
print(f"\nReplicate split-half consistency: Pearson r = {r_rep:+.4f} (p = {p_rep:.4f})")

# Correlation with mean response (C0, C1, C2)
r_c0, p_c0 = pearsonr(df_aissa['arti'], df_aissa['c0_pseudobulk_l2'])
r_c1, p_c1 = pearsonr(df_aissa['arti'], df_aissa['c1_energy_dist'])
r_c2, p_c2 = pearsonr(df_aissa['arti'], df_aissa['c2_ps_magnitude'])
print(f"Correlation r(ARTI, C0 Pseudobulk L2) = {r_c0:+.4f} (p = {p_c0:.4e})")
print(f"Correlation r(ARTI, C1 Energy Dist)   = {r_c1:+.4f} (p = {p_c1:.4e})")
print(f"Correlation r(ARTI, C2 PS Magnitude)  = {r_c2:+.4f} (p = {p_c2:.4e})")

t1 = time.time()
print(f"Completed Aissa et al. benchmark in {t1-t0:.1f}s.")
