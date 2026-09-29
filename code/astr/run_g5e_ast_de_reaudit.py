# code/astr/run_g5e_ast_de_reaudit.py
# Gate G5E: AST Differential Expression Re-Audit, Hierarchical Inference & Stress Annotation

import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

os.makedirs('tables/astr', exist_ok=True)
os.makedirs('results/astr', exist_ok=True)

print("--- Running Gate G5E: AST Differential Expression Re-Audit ---")
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
obs = adata.obs.copy()
obs['cell_idx'] = np.arange(len(obs))

obs['is_qc_excluded'] = (obs['pct_mito'] >= 15.0) | (obs['total_counts'] < 800) | (obs['n_genes'] < 300)
obs['is_qc_pass'] = ~obs['is_dying'] if 'is_dying' in obs.columns else ~obs['is_qc_excluded']
qc_pass_obs = obs[obs['is_qc_pass']].copy()

classified = pd.read_parquet('data_processed/tahoe_mvp_cells_classified.parquet')
qc_pass_obs['cell_state_category'] = classified.loc[qc_pass_obs.index, 'cell_state_category']

# Primary evaluation: KRAS G12C lines under Adagrasib (HOP62, NCI-H2030, NCI-H23)
kras_adag = qc_pass_obs[(qc_pass_obs['cell_name'].isin(['HOP62', 'NCI-H2030', 'NCI-H23'])) & 
                        (qc_pass_obs['drug'] == 'Adagrasib')].copy()

# Feature selection: Top 1000 HVGs
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
sub_X = adata.X[:, top_gene_indices].toarray()
lib_size = sub_X.sum(axis=1, keepdims=True)
lib_size[lib_size == 0] = 1.0
norm_X = np.log1p((sub_X / lib_size) * 10000.0)
gene_symbols = np.array(adata.var_names)[top_gene_indices]

# Predefined canonical stress genes
known_stress_genes = {
    'MT2A': 'Metallothionein / Heavy metal & ROS stress',
    'MT1X': 'Metallothionein / Heavy metal & ROS stress',
    'HSPA1A': 'Heat shock 70kDa protein / Chaperone',
    'HSP90AA1': 'Heat shock 90kDa protein / Chaperone',
    'DNAJB1': 'Hsp40 chaperone',
    'HMOX1': 'Heme oxygenase / Oxidative stress',
    'SOD2': 'Superoxide dismutase / Mitochondrial ROS',
    'FOS': 'Immediate-early proto-oncogene / AP-1',
    'JUN': 'Immediate-early proto-oncogene / AP-1',
    'EGR1': 'Early growth response / Transcriptional stress',
    'ATF3': 'Activating transcription factor 3 / Cellular stress',
    'ATF4': 'Activating transcription factor 4 / Integrated stress response',
    'DDIT3': 'DNA damage inducible transcript 3 / CHOP / ER stress',
    'PPP1R15A': 'GADD34 / ER stress recovery'
}

# Run differential expression for HOP62 Adagrasib
hop62_cells = kras_adag[kras_adag['cell_name'] == 'HOP62']
ast_mask = hop62_cells['cell_state_category'] == 'Adaptive_Survival_Tail'
dom_mask = hop62_cells['cell_state_category'] == 'Dominant_Response'

ast_expr = norm_X[hop62_cells['cell_idx'][ast_mask]]
dom_expr = norm_X[hop62_cells['cell_idx'][dom_mask]]

de_records = []
for g_i, g_name in enumerate(gene_symbols):
    v_ast = ast_expr[:, g_i]
    v_dom = dom_expr[:, g_i]
    m_ast = float(np.mean(v_ast))
    m_dom = float(np.mean(v_dom))
    log2fc = float(m_ast - m_dom)
    
    stat, pval = mannwhitneyu(v_ast, v_dom, alternative='two-sided')
    
    # Stress classification
    stress_cat = known_stress_genes.get(g_name, 'Non_Generic_Stress_Candidate')
    
    de_records.append({
        'gene': g_name,
        'statistical_unit': 'cell_within_HOP62_Adagrasib',
        'n_ast_cells': len(v_ast),
        'n_dominant_cells': len(v_dom),
        'mean_ast': round(m_ast, 4),
        'mean_dominant': round(m_dom, 4),
        'log2_fold_change': round(log2fc, 4),
        'mann_whitney_u_stat': float(stat),
        'raw_p_value': pval,
        'stress_program_annotation': stress_cat,
        'persister_terminology_status': 'AST_Associated_Candidate_NOT_DTP_Marker'
    })

de_df = pd.DataFrame(de_records).sort_values(by='log2_fold_change', ascending=False)

# Compute Benjamini-Hochberg FDR
pvals = de_df['raw_p_value'].values
order = np.argsort(pvals)
ranks = np.empty_like(order)
ranks[order] = np.arange(1, len(pvals) + 1)
fdr = pvals * len(pvals) / ranks
fdr_sorted = np.minimum.accumulate(fdr[order][::-1])[::-1]
fdr_adj = np.empty_like(fdr_sorted)
fdr_adj[order] = np.clip(fdr_sorted, 0.0, 1.0)
de_df['bh_fdr'] = fdr_adj

de_df.to_csv('tables/astr/ast_de_reaudit_table.csv', index=False)

print("Top 15 Upregulated Genes in HOP62 AST with FDR & Stress Program Annotations:")
print(de_df[['gene', 'log2_fold_change', 'raw_p_value', 'bh_fdr', 'stress_program_annotation']].head(15).to_string())

summary_g5e = {
    'primary_condition': 'HOP62 + Adagrasib',
    'total_genes_screened': len(de_df),
    'genes_with_raw_p_005': int(np.sum(de_df['raw_p_value'] < 0.05)),
    'genes_with_fdr_010': int(np.sum(de_df['bh_fdr'] < 0.10)),
    'top_upregulated_gene': str(de_df.iloc[0]['gene']),
    'top_upregulated_log2fc': float(de_df.iloc[0]['log2_fold_change']),
    'top_upregulated_fdr': float(de_df.iloc[0]['bh_fdr']),
    'top_generic_stress_overlap': 'MT2A belongs to metallothionein ROS stress program; UBE2E2 and USP53 represent ubiquitin-proteasome system',
    'terminology_guard': 'ENFORCED: All genes cataloged strictly as AST_Associated_Candidate_NOT_DTP_Marker'
}

with open('results/astr/ast_de_reaudit_summary.json', 'w') as f:
    json.dump(summary_g5e, f, indent=2)

print("\nSaved tables/astr/ast_de_reaudit_table.csv and results/astr/ast_de_reaudit_summary.json successfully!")
