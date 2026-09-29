import os
import json
import anndata
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

print('--- Testing Pan-Condition DTP Similarity across all High-AST Perturbations ---')
tahoe_adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
classified_df = pd.read_parquet('data_processed/tahoe_mvp_cells_classified.parquet')
tahoe_adata.obs = classified_df

with open('results/external_dtp/dtp_reference_signatures.json') as f:
    dtp_ref = json.load(f)

dtp_de_df = pd.read_csv('tables/external_dtp/dtp_signature_genes.csv')
gse_adata = anndata.read_h5ad('data_processed/external_dtp/gse275199_processed_adata.h5ad')

shared_genes = sorted(list(set(tahoe_adata.var_names).intersection(set(gse_adata.var_names))))
dtp_de_shared = dtp_de_df[dtp_de_df['gene'].isin(shared_genes)].copy()
top_hvg_genes = dtp_de_shared.sort_values(by='log2FC', key=abs, ascending=False).head(1000)['gene'].tolist()
dtp_vector = dtp_de_shared.set_index('gene').loc[top_hvg_genes, 'log2FC'].values

up_genes = [g for g in dtp_ref['top100_dtp_upregulated'] if g in tahoe_adata.var_names]
down_genes = [g for g in dtp_ref['top100_dtp_downregulated'] if g in tahoe_adata.var_names]

# Test top conditions with AST >= 5 cells
ast_counts = tahoe_adata.obs[tahoe_adata.obs['cell_state_category'] == 'Adaptive_Survival_Tail'].groupby(['cell_name', 'drug', 'dose', 'plate'], observed=True).size()
top_ast_conds = ast_counts[ast_counts >= 4].reset_index()

pan_results = []
for idx, r in top_ast_conds.iterrows():
    cname, drug, dose, plate = r['cell_name'], r['drug'], r['dose'], r['plate']
    cond_cells = tahoe_adata[(tahoe_adata.obs['cell_name'] == cname) & 
                             (tahoe_adata.obs['drug'] == drug) & 
                             (tahoe_adata.obs['dose'] == dose) & 
                             (tahoe_adata.obs['plate'] == plate) & 
                             (tahoe_adata.obs['is_viable'])].copy()
    
    dmso_cells = tahoe_adata[(tahoe_adata.obs['cell_name'] == cname) & 
                             (tahoe_adata.obs['drug'] == 'DMSO_TF') & 
                             (tahoe_adata.obs['plate'] == plate) & 
                             (tahoe_adata.obs['is_viable'])].copy()
    
    if len(dmso_cells) < 10:
        continue
        
    def get_norm(sub_ad, genes):
        sub = sub_ad[:, genes].X
        if hasattr(sub, 'toarray'):
            sub = sub.toarray()
        lib = sub_ad.obs['total_counts'].values[:, None]
        lib[lib == 0] = 1.0
        return np.log1p((sub / lib) * 1e4)
        
    cond_norm = get_norm(cond_cells, top_hvg_genes)
    dmso_norm = get_norm(dmso_cells, top_hvg_genes)
    dmso_mean = np.mean(dmso_norm, axis=0)
    delta_expr = cond_norm - dmso_mean
    
    rhos = [spearmanr(delta_expr[i], dtp_vector)[0] for i in range(len(delta_expr))]
    
    up_m = get_norm(cond_cells, up_genes)
    down_m = get_norm(cond_cells, down_genes)
    dtp_scores = np.mean(up_m, axis=1) - np.mean(down_m, axis=1)
    
    ast_mask = (cond_cells.obs['cell_state_category'] == 'Adaptive_Survival_Tail').values
    dom_mask = (cond_cells.obs['cell_state_category'] == 'Dominant_Response').values
    
    if np.sum(ast_mask) < 3 or np.sum(dom_mask) < 10:
        continue
        
    ast_rhos = np.array(rhos)[ast_mask]
    dom_rhos = np.array(rhos)[dom_mask]
    
    ast_scores = np.array(dtp_scores)[ast_mask]
    dom_scores = np.array(dtp_scores)[dom_mask]
    
    diff_rho = np.mean(ast_rhos) - np.mean(dom_rhos)
    u_rho, p_rho = mannwhitneyu(ast_rhos, dom_rhos, alternative='greater')
    
    diff_score = np.mean(ast_scores) - np.mean(dom_scores)
    u_score, p_score = mannwhitneyu(ast_scores, dom_scores, alternative='greater')
    
    pan_results.append({
        'cell_name': cname,
        'drug': drug,
        'dose': dose,
        'n_ast': int(np.sum(ast_mask)),
        'n_dom': int(np.sum(dom_mask)),
        'ast_mean_rho': round(float(np.mean(ast_rhos)), 4),
        'dom_mean_rho': round(float(np.mean(dom_rhos)), 4),
        'delta_rho': round(float(diff_rho), 4),
        'p_val_rho': round(float(p_rho), 4),
        'ast_mean_score': round(float(np.mean(ast_scores)), 4),
        'dom_mean_score': round(float(np.mean(dom_scores)), 4),
        'delta_score': round(float(diff_score), 4),
        'p_val_score': round(float(p_score), 4)
    })

pan_df = pd.DataFrame(pan_results)
print('\n--- Pan-Condition DTP Projection Summary ---')
print(pan_df.to_string())
pan_df.to_csv('tables/external_dtp/pan_condition_dtp_projections.csv', index=False)
