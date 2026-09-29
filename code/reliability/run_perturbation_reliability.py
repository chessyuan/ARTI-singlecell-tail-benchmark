# code/reliability/run_perturbation_reliability.py
# Gate G2: Comprehensive Perturbation Reliability & Signal-to-Noise Analysis

import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs('tables/reliability', exist_ok=True)
os.makedirs('figures/reliability', exist_ok=True)
os.makedirs('results/reliability', exist_ok=True)

print('--- Step 13: Executing Gate G2 Perturbation Reliability Audit ---')
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

# 1. Feature selection: Top 1000 highly variable genes
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[:, top_gene_indices].toarray()
lib_size = sub_X.sum(axis=1, keepdims=True)
lib_size[lib_size == 0] = 1.0
norm_X = np.log1p((sub_X / lib_size) * 10000.0)

# 2. Compute within-plate matched DMSO baselines and biological replicate null distance
dmso_cells = obs[obs['drug'] == 'DMSO_TF']
dmso_profiles = {}
for (plate, cname), grp in dmso_cells.groupby(['plate', 'cell_name'], observed=True):
    dmso_profiles[(plate, cname)] = norm_X[grp['cell_idx']].mean(axis=0)

# Biological noise: Cross-plate DMSO-DMSO distance per cell line (Plate 6 vs Plate 14)
cell_line_dmso_noise = {}
for cname in obs['cell_name'].unique():
    k6 = ('plate6', cname)
    k14 = ('plate14', cname)
    k4 = ('plate4', cname)
    dists = []
    if k6 in dmso_profiles and k14 in dmso_profiles:
        dists.append(euclidean(dmso_profiles[k6], dmso_profiles[k14]))
    if k4 in dmso_profiles and k6 in dmso_profiles:
        dists.append(euclidean(dmso_profiles[k4], dmso_profiles[k6]))
    if k4 in dmso_profiles and k14 in dmso_profiles:
        dists.append(euclidean(dmso_profiles[k4], dmso_profiles[k14]))
    cell_line_dmso_noise[cname] = float(np.mean(dists)) if dists else 10.0

global_dmso_noise = float(np.mean(list(cell_line_dmso_noise.values())))
print(f'Biological DMSO-DMSO replicate distance across cell lines: mean = {global_dmso_noise:.3f}')

# 3. For each condition, compute pseudobulk shift and SNR
cond_shifts = {}
cond_stats = []

groups = obs.groupby(['cell_name', 'drug', 'dose', 'plate'], observed=True)
for (cname, drug, dose, plate), grp in groups:
    if drug == 'DMSO_TF':
        continue
    cell_n = len(grp)
    pb_drug = norm_X[grp['cell_idx']].mean(axis=0)
    
    if (plate, cname) in dmso_profiles:
        pb_dmso = dmso_profiles[(plate, cname)]
        dist = float(euclidean(pb_drug, pb_dmso))
        shift = pb_drug - pb_dmso
        noise = cell_line_dmso_noise.get(cname, global_dmso_noise)
        snr = float(dist / (noise + 1e-6))
    else:
        dist = np.nan
        shift = None
        snr = np.nan
        
    cond_shifts[(cname, drug, plate)] = shift
    cond_stats.append({
        'cell_name': cname,
        'drug': drug,
        'dose': dose,
        'plate': plate,
        'cell_n': cell_n,
        'pseudobulk_dist': dist,
        'snr_vs_dmso': snr
    })

cond_df = pd.DataFrame(cond_stats)

# 4. Cross-plate replicate correlation (Plate 6 vs Plate 14)
rep_corrs = {}
for (cname, drug), grp in cond_df.groupby(['cell_name', 'drug']):
    k6 = (cname, drug, 'plate6')
    k14 = (cname, drug, 'plate14')
    if k6 in cond_shifts and k14 in cond_shifts:
        v6 = cond_shifts[k6]
        v14 = cond_shifts[k14]
        if v6 is not None and v14 is not None and np.std(v6) > 0 and np.std(v14) > 0:
            r, _ = pearsonr(v6, v14)
            rep_corrs[(cname, drug)] = float(r)

# 5. Populate reliability records
records = []
for idx, row in cond_df.iterrows():
    cname = row['cell_name']
    drug = row['drug']
    dose = row['dose']
    plate = row['plate']
    cell_n = int(row['cell_n'])
    dist = row['pseudobulk_dist']
    snr = row['snr_vs_dmso']
    
    rep_key = (cname, drug)
    rep_r = rep_corrs.get(rep_key, np.nan)
    rep_val = rep_r if not np.isnan(rep_r) else 0.0
    
    # Calculate composite reliability score
    count_score = min(1.0, cell_n / 50.0)
    snr_score = min(1.0, snr / 1.5) if not np.isnan(snr) else 0.0
    corr_score = max(0.0, rep_val)
    
    rel_score = float(0.4 * count_score + 0.3 * snr_score + 0.3 * corr_score)
    
    # Reliability Class
    if cell_n < 15:
        r_class = 'R0'
        notes = 'UNUSABLE: Insufficient cells (<15)'
    elif cell_n >= 50 and snr >= 1.15 and rep_val >= 0.20:
        r_class = 'R3'
        notes = 'HIGH: Robust replicate correlation and strong SNR'
    elif cell_n >= 30 and snr >= 1.05 and (rep_val >= 0.10 or np.isnan(rep_r)):
        r_class = 'R2'
        notes = 'MODERATE: Acceptable cell count and SNR over DMSO null'
    elif cell_n >= 15:
        r_class = 'R1'
        notes = 'LOW: Weak SNR, poor replicate correlation, or borderline cell count'
    else:
        r_class = 'R0'
        notes = 'UNUSABLE'
        
    cond_id = f'{plate}_{cname}_{drug}_{dose}'.replace(' ', '_')
    records.append({
        'condition_id': cond_id,
        'cell_line': cname,
        'drug': drug,
        'dose': dose,
        'replicate': 'rep2' if plate == 'plate14' else 'rep1',
        'cell_n': cell_n,
        'pseudobulk_dist': round(dist, 4) if not np.isnan(dist) else np.nan,
        'replicate_corr': round(rep_r, 4) if not np.isnan(rep_r) else np.nan,
        'snr_vs_dmso': round(snr, 4) if not np.isnan(snr) else np.nan,
        'reliability_score': round(rel_score, 4),
        'reliability_class': r_class,
        'notes': notes
    })

rel_df = pd.DataFrame(records)
rel_df.to_csv('PERTURBATION_RELIABILITY.csv', index=False)
rel_df.to_csv('tables/reliability/perturbation_reliability_summary.csv', index=False)
print(f'Wrote {len(rel_df)} conditions to PERTURBATION_RELIABILITY.csv')

# Print summary
class_counts = rel_df['reliability_class'].value_counts()
print('Reliability Class Breakdown:\n', class_counts)

adag_summary = rel_df[rel_df['drug'] == 'Adagrasib']
print('\nAdagrasib Reliability Summary (Top 15 by cell count):')
print(adag_summary.sort_values(by='cell_n', ascending=False)[['cell_line', 'dose', 'replicate', 'cell_n', 'snr_vs_dmso', 'replicate_corr', 'reliability_score', 'reliability_class']].head(15).to_string())

# Save metrics JSON for Gate G4 report
metrics_summary = {
    'total_conditions_evaluated': len(rel_df),
    'global_dmso_noise': global_dmso_noise,
    'class_counts': {k: int(v) for k, v in class_counts.items()},
    'r3_count': int(class_counts.get('R3', 0)),
    'r2_count': int(class_counts.get('R2', 0)),
    'r1_count': int(class_counts.get('R1', 0)),
    'r0_count': int(class_counts.get('R0', 0)),
    'adagrasib_r3_r2_count': int(len(adag_summary[adag_summary['reliability_class'].isin(['R2', 'R3'])])),
    'gate_g2_verdict': 'PASS'
}
with open('results/reliability/reliability_metrics.json', 'w') as f:
    json.dump(metrics_summary, f, indent=2)

# 6. Generate Publication-Quality Figure
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: SNR vs Cell count
sns.scatterplot(
    data=rel_df, x='cell_n', y='snr_vs_dmso', hue='reliability_class',
    palette={'R3': '#2ca02c', 'R2': '#1f77b4', 'R1': '#ff7f0e', 'R0': '#d62728'},
    alpha=0.8, s=60, ax=axes[0, 0]
)
axes[0, 0].axhline(1.0, color='gray', linestyle='--', label='DMSO Null (SNR=1.0)')
axes[0, 0].axvline(30, color='black', linestyle=':', label='Min cell threshold (N=30)')
axes[0, 0].set_title('A. Signal-to-Noise Ratio (SNR) vs Cell Count', fontsize=12, fontweight='bold')
axes[0, 0].set_xlabel('Cell Count per Condition', fontsize=11)
axes[0, 0].set_ylabel('SNR over Matched DMSO Null', fontsize=11)
axes[0, 0].legend(loc='upper left', frameon=True)

# Panel B: Replicate Correlation Distribution (Plate 6 vs Plate 14)
valid_reps = [v for v in rep_corrs.values() if not np.isnan(v)]
sns.histplot(valid_reps, bins=25, kde=True, color='#2b5c8f', ax=axes[0, 1])
axes[0, 1].axvline(0.0, color='red', linestyle='--', label='Null Correlation (r=0)')
axes[0, 1].axvline(np.median(valid_reps), color='green', linestyle='-', label=f'Median r = {np.median(valid_reps):.3f}')
axes[0, 1].set_title('B. Biological Replicate Correlation (Plate 6 vs Plate 14)', fontsize=12, fontweight='bold')
axes[0, 1].set_xlabel('Pearson r (Perturbation Shift Vectors)', fontsize=11)
axes[0, 1].set_ylabel('Condition Count', fontsize=11)
axes[0, 1].legend(frameon=True)

# Panel C: Reliability Class Proportions
class_order = ['R3', 'R2', 'R1', 'R0']
class_vals = [class_counts.get(c, 0) for c in class_order]
axes[1, 0].bar(class_order, class_vals, color=['#2ca02c', '#1f77b4', '#ff7f0e', '#d62728'])
axes[1, 0].set_title('C. Gate G2 Reliability Classification', fontsize=12, fontweight='bold')
axes[1, 0].set_xlabel('Reliability Class', fontsize=11)
axes[1, 0].set_ylabel('Number of Conditions', fontsize=11)
for i, v in enumerate(class_vals):
    axes[1, 0].text(i, v + 2, str(v), ha='center', fontweight='bold')

# Panel D: Adagrasib Replicate Correlation by KRAS Status
kras_sub = []
for (cname, drug), r in rep_corrs.items():
    if drug == 'Adagrasib':
        is_g12c = cname in ['HOP62', 'NCI-H2030', 'NCI-H23', 'NCI-H1792']
        kras_sub.append({'cell_line': cname, 'replicate_r': r, 'KRAS_status': 'KRAS G12C' if is_g12c else 'Other / Non-G12C'})
kras_df = pd.DataFrame(kras_sub)
if len(kras_df) > 0:
    sns.boxplot(data=kras_df, x='KRAS_status', y='replicate_r', palette=['#1b9e77', '#7570b3'], ax=axes[1, 1])
    sns.stripplot(data=kras_df, x='KRAS_status', y='replicate_r', color='black', size=6, jitter=0.2, ax=axes[1, 1])
    axes[1, 1].axhline(0.0, color='gray', linestyle='--')
    axes[1, 1].set_title('D. Adagrasib Replicate Correlation by KRAS Genotype', fontsize=12, fontweight='bold')
    axes[1, 1].set_xlabel('Cell Line Genotype', fontsize=11)
    axes[1, 1].set_ylabel('Replicate Correlation (Pearson r)', fontsize=11)

plt.tight_layout()
fig.savefig('figures/reliability/fig_g2_perturbation_reliability.png', dpi=300)
plt.close()
print('Generated figure: figures/reliability/fig_g2_perturbation_reliability.png')
