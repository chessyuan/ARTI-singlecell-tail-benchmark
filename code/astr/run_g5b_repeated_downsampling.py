# code/astr/run_g5b_repeated_downsampling.py
# Gate G5B: Repeated Downsampling Stability across Multiple Sample Sizes & Seeds

import os
import json
import anndata
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import pearsonr, spearmanr

os.makedirs('tables/astr', exist_ok=True)
os.makedirs('results/astr', exist_ok=True)

print("--- Running Gate G5B: Repeated Downsampling Stability Protocol ---")
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

obs['is_qc_excluded'] = (obs['pct_mito'] >= 15.0) | (obs['total_counts'] < 800) | (obs['n_genes'] < 300)
obs['is_qc_pass'] = ~obs['is_qc_excluded']
qc_pass_obs = obs[obs['is_qc_pass']].copy()

# Top 1000 HVGs
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[:, top_gene_indices].toarray()
lib_size = sub_X.sum(axis=1, keepdims=True)
lib_size[lib_size == 0] = 1.0
norm_X = np.log1p((sub_X / lib_size) * 10000.0)

# Full condition results
astr_full = pd.read_csv('tables/astr/astr_all_conditions.csv')
dmso_map = dict(zip(zip(astr_full['cell_name'], astr_full['plate']), astr_full['dmso_baseline_dispersion']))

# Filter eligible conditions: N >= 50
eligible_30 = astr_full[astr_full['n_qc_pass'] >= 50].copy()
eligible_50 = astr_full[astr_full['n_qc_pass'] >= 80].copy()

print(f"Eligible conditions for N=30 downsampling (N_qc >= 50): {len(eligible_30)}")
print(f"Eligible conditions for N=50 downsampling (N_qc >= 80): {len(eligible_50)}")

seeds = [42 + i * 17 for i in range(50)]

results_by_n = {}
for target_n, sub_df, label in [(30, eligible_30, 'N30'), (50, eligible_50, 'N50')]:
    if len(sub_df) < 10:
        continue
    
    pearson_list = []
    spearman_list = []
    top_decile_overlap_list = []
    class_agreement_list = []
    
    n_top = max(3, int(len(sub_df) * 0.10))
    full_top_conditions = set(sub_df.sort_values(by='delta_dispersion', ascending=False)['condition_id'].iloc[:n_top])
    
    for s_idx, seed in enumerate(seeds):
        np.random.seed(seed)
        ds_deltas = []
        cond_ids = []
        
        for _, r in sub_df.iterrows():
            cid = r['condition_id']
            cname = r['cell_name']
            drug = r['drug']
            plate = r['plate']
            base_disp = dmso_map.get((cname, plate), np.nan)
            
            c_cells = qc_pass_obs[(qc_pass_obs['cell_name'] == cname) & (qc_pass_obs['drug'] == drug) & (qc_pass_obs['plate'] == plate)]
            sampled_idx = np.random.choice(c_cells['cell_idx'].values, size=target_n, replace=False)
            
            pca = PCA(n_components=min(5, target_n - 1), random_state=seed)
            pcs = pca.fit_transform(norm_X[sampled_idx])
            centroid = np.median(pcs, axis=0, keepdims=True)
            dists = np.linalg.norm(pcs - centroid, axis=1)
            ds_disp = float(np.mean(dists))
            delta = ds_disp - base_disp
            
            ds_deltas.append(delta)
            cond_ids.append(cid)
            
        r_p, _ = pearsonr(sub_df['delta_dispersion'], ds_deltas)
        r_s, _ = spearmanr(sub_df['delta_dispersion'], ds_deltas)
        pearson_list.append(r_p)
        spearman_list.append(r_s)
        
        # Top decile overlap
        res_temp = pd.DataFrame({'condition_id': cond_ids, 'delta': ds_deltas})
        ds_top = set(res_temp.sort_values(by='delta', ascending=False)['condition_id'].iloc[:n_top])
        overlap = len(full_top_conditions.intersection(ds_top)) / n_top
        top_decile_overlap_list.append(overlap)
        
        # Agreement on operational threshold (delta > 0.50)
        full_sig = sub_df['delta_dispersion'] > 0.50
        ds_sig = np.array(ds_deltas) > 0.50
        agreement = float(np.mean(full_sig.values == ds_sig))
        class_agreement_list.append(agreement)
        
    results_by_n[label] = {
        'target_n': target_n,
        'n_conditions': len(sub_df),
        'n_seeds': len(seeds),
        'pearson_median': round(float(np.median(pearson_list)), 4),
        'pearson_iqr': round(float(np.percentile(pearson_list, 75) - np.percentile(pearson_list, 25)), 4),
        'pearson_ci95': [round(float(np.percentile(pearson_list, 2.5)), 4), round(float(np.percentile(pearson_list, 97.5)), 4)],
        'spearman_median': round(float(np.median(spearman_list)), 4),
        'spearman_iqr': round(float(np.percentile(spearman_list, 75) - np.percentile(spearman_list, 25)), 4),
        'spearman_ci95': [round(float(np.percentile(spearman_list, 2.5)), 4), round(float(np.percentile(spearman_list, 97.5)), 4)],
        'top_decile_overlap_median': round(float(np.median(top_decile_overlap_list)), 4),
        'top_decile_overlap_ci95': [round(float(np.percentile(top_decile_overlap_list, 2.5)), 4), round(float(np.percentile(top_decile_overlap_list, 97.5)), 4)],
        'class_agreement_median': round(float(np.median(class_agreement_list)), 4),
        'class_agreement_ci95': [round(float(np.percentile(class_agreement_list, 2.5)), 4), round(float(np.percentile(class_agreement_list, 97.5)), 4)],
        'raw_pearson': pearson_list,
        'raw_spearman': spearman_list
    }

# Save table of seed runs
stability_rows = []
for label, data in results_by_n.items():
    stability_rows.append({
        'subset_size': data['target_n'],
        'n_conditions': data['n_conditions'],
        'pearson_median': data['pearson_median'],
        'pearson_iqr': data['pearson_iqr'],
        'pearson_ci95_low': data['pearson_ci95'][0],
        'pearson_ci95_high': data['pearson_ci95'][1],
        'spearman_median': data['spearman_median'],
        'spearman_iqr': data['spearman_iqr'],
        'spearman_ci95_low': data['spearman_ci95'][0],
        'spearman_ci95_high': data['spearman_ci95'][1],
        'top_decile_overlap_median': data['top_decile_overlap_median'],
        'top_decile_overlap_ci95_low': data['top_decile_overlap_ci95'][0],
        'top_decile_overlap_ci95_high': data['top_decile_overlap_ci95'][1],
        'class_agreement_median': data['class_agreement_median'],
        'class_agreement_ci95_low': data['class_agreement_ci95'][0],
        'class_agreement_ci95_high': data['class_agreement_ci95'][1]
    })

stab_df = pd.DataFrame(stability_rows)
stab_df.to_csv('tables/astr/ast_downsampling_repeated_stability.csv', index=False)

summary_clean = {k: {sub_k: v for sub_k, v in data.items() if not sub_k.startswith('raw_')} for k, data in results_by_n.items()}
with open('results/astr/downsampling_stability_summary.json', 'w') as f:
    json.dump(summary_clean, f, indent=2)

print("--- Gate G5B Repeated Downsampling Summary ---")
print(json.dumps(summary_clean, indent=2))
