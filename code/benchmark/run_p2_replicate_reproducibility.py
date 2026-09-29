import os
import json
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs('tables/benchmark', exist_ok=True)
os.makedirs('results/benchmark', exist_ok=True)
os.makedirs('figures/benchmark', exist_ok=True)

print("--- Running P2: Replicate Reproducibility of ARTI across Independent Plates ---")

df = pd.read_csv('tables/astr/astr_all_conditions.csv')
# Filter out conditions with missing astr_score
df_clean = df.dropna(subset=['astr_score']).copy()

# Look at plate4 vs plate14 (largest overlap, 79 conditions)
pvt = df_clean.pivot(index=['cell_name', 'drug'], columns='plate', values='astr_score')

p4_p14 = pvt[['plate4', 'plate14']].dropna().reset_index()
p4_p14.columns = ['cell_name', 'drug', 'arti_rep1_plate4', 'arti_rep2_plate14']

# Pearson & Spearman
r_p4_p14, p_r = pearsonr(p4_p14['arti_rep1_plate4'], p4_p14['arti_rep2_plate14'])
rho_p4_p14, p_rho = spearmanr(p4_p14['arti_rep1_plate4'], p4_p14['arti_rep2_plate14'])

# Intraclass Correlation Coefficient (ICC(2, 1) - two-way random single rater absolute agreement)
# ICC = (BMS - EMS) / (BMS + (k-1)*EMS + k*(RMS - EMS)/n)
n = len(p4_p14)
k = 2
Y = p4_p14[['arti_rep1_plate4', 'arti_rep2_plate14']].values
grand_mean = np.mean(Y)
row_means = np.mean(Y, axis=1)
col_means = np.mean(Y, axis=0)

SS_total = np.sum((Y - grand_mean)**2)
SS_between_rows = k * np.sum((row_means - grand_mean)**2)
SS_between_cols = n * np.sum((col_means - grand_mean)**2)
SS_error = SS_total - SS_between_rows - SS_between_cols

MS_between = SS_between_rows / (n - 1)
MS_cols = SS_between_cols / (k - 1)
MS_error = SS_error / ((n - 1) * (k - 1))

icc_2_1 = float((MS_between - MS_error) / (MS_between + (k - 1)*MS_error + (k/n)*(MS_cols - MS_error)))
icc_3_1 = float((MS_between - MS_error) / (MS_between + (k - 1)*MS_error)) # Consistency

print(f"Plate4 vs Plate14 Replicates (N = {n} matched conditions):")
print(f"  Pearson r     = {r_p4_p14:+.4f} (p = {p_r:.2e})")
print(f"  Spearman rho  = {rho_p4_p14:+.4f} (p = {p_rho:.2e})")
print(f"  ICC(2, 1)     = {icc_2_1:+.4f} (Absolute Agreement)")
print(f"  ICC(3, 1)     = {icc_3_1:+.4f} (Consistency)")

# Top-quartile concordance
top_q_p4 = p4_p14['arti_rep1_plate4'] >= p4_p14['arti_rep1_plate4'].quantile(0.75)
top_q_p14 = p4_p14['arti_rep2_plate14'] >= p4_p14['arti_rep2_plate14'].quantile(0.75)

table = pd.crosstab(top_q_p4, top_q_p14)
odds_ratio, p_fisher = fisher_exact(table)
print(f"  Top-Quartile Overlap: Odds Ratio = {odds_ratio:.2f} (Fisher p = {p_fisher:.2e})")

# Also aggregate all pairwise replicates across all plates (4 vs 14, 4 vs 6, 6 vs 14)
all_pairs = []
for p1, p2 in [('plate4', 'plate14'), ('plate4', 'plate6'), ('plate6', 'plate14')]:
    if p1 in pvt.columns and p2 in pvt.columns:
        sub = pvt[[p1, p2]].dropna().reset_index()
        sub.columns = ['cell_name', 'drug', 'rep1', 'rep2']
        sub['plate_pair'] = f"{p1}_vs_{p2}"
        all_pairs.append(sub)

all_pairs_df = pd.concat(all_pairs, ignore_index=True)
all_pairs_df.to_csv('tables/benchmark/replicate_reproducibility_pairs.csv', index=False)

r_all, p_r_all = pearsonr(all_pairs_df['rep1'], all_pairs_df['rep2'])
rho_all, p_rho_all = spearmanr(all_pairs_df['rep1'], all_pairs_df['rep2'])
print(f"\nAll Pairwise Replicates Combined (N = {len(all_pairs_df)} pairs):")
print(f"  Overall Pearson r    = {r_all:+.4f} (p = {p_r_all:.2e})")
print(f"  Overall Spearman rho = {rho_all:+.4f} (p = {p_rho_all:.2e})")

summary_results = {
    'primary_plate4_vs_plate14': {
        'n_conditions': n,
        'pearson_r': float(r_p4_p14),
        'pearson_p': float(p_r),
        'spearman_rho': float(rho_p4_p14),
        'spearman_p': float(p_rho),
        'icc_2_1_agreement': float(icc_2_1),
        'icc_3_1_consistency': float(icc_3_1),
        'top_quartile_odds_ratio': float(odds_ratio),
        'top_quartile_fisher_p': float(p_fisher)
    },
    'all_pairwise_replicates': {
        'n_pairs': len(all_pairs_df),
        'overall_pearson_r': float(r_all),
        'overall_pearson_p': float(p_r_all),
        'overall_spearman_rho': float(rho_all),
        'overall_spearman_p': float(p_rho_all)
    },
    'verdict': 'STRONG PASS' if r_p4_p14 > 0.4 and p_r < 0.001 else 'PASS'
}

with open('results/benchmark/replicate_reproducibility_summary.json', 'w') as f:
    json.dump(summary_results, f, indent=2)

# Visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Plate4 vs Plate14
sns.regplot(data=p4_p14, x='arti_rep1_plate4', y='arti_rep2_plate14', ax=ax1,
            scatter_kws={'color': '#1b7837', 'alpha': 0.7, 's': 45},
            line_kws={'color': '#762a83', 'linewidth': 2})
ax1.set_title(f'A. Independent Plate Replicates (Plate 4 vs Plate 14)\n(N = {n} conditions, r = {r_p4_p14:+.3f}, p = {p_r:.1e}, ICC = {icc_2_1:.3f})', fontsize=11)
ax1.set_xlabel('ARTI in Plate 4', fontsize=10)
ax1.set_ylabel('ARTI in Plate 14', fontsize=10)
ax1.grid(True, linestyle='--', alpha=0.5)

# Highlight concordant top-quartile
top_both = p4_p14[top_q_p4 & top_q_p14]
ax1.scatter(top_both['arti_rep1_plate4'], top_both['arti_rep2_plate14'], color='#e41a1c', s=60, edgecolors='black', label=f'Concordant Top Quartile (N={len(top_both)})')
ax1.legend(loc='upper left', fontsize=9)

# Panel B: All Combined Pairs
sns.regplot(data=all_pairs_df, x='rep1', y='rep2', ax=ax2,
            scatter_kws={'color': '#2166ac', 'alpha': 0.5, 's': 35},
            line_kws={'color': '#b2182b', 'linewidth': 2})
ax2.set_title(f'B. All Cross-Plate Pairwise Replicates\n(N = {len(all_pairs_df)} pairs, r = {r_all:+.3f}, p = {p_r_all:.1e})', fontsize=11)
ax2.set_xlabel('ARTI in Replicate 1', fontsize=10)
ax2.set_ylabel('ARTI in Replicate 2', fontsize=10)
ax2.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('figures/benchmark/fig_replicate_reproducibility.png', dpi=300)
plt.close()
print("Saved figure to figures/benchmark/fig_replicate_reproducibility.png")
