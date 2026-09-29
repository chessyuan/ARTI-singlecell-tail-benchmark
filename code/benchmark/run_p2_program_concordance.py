import anndata
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

print("--- Testing State-Program Concordance across Independent Replicates ---")
adata = anndata.read_h5ad('data_processed/tahoe_mvp_adata.h5ad')
classified = pd.read_parquet('data_processed/tahoe_mvp_cells_classified.parquet')

# Join classification to adata
adata.obs = adata.obs.set_index('cell_id')
adata.obs['cell_state_category'] = classified.set_index('cell_id')['cell_state_category']

# Find conditions with ART cells in both Plate 4 and Plate 14
p4_p14_matched = pd.read_csv('tables/benchmark/replicate_reproducibility_pairs.csv')
p4_p14_matched = p4_p14_matched[p4_p14_matched['plate_pair'] == 'plate4_vs_plate14']

# 1000 HVGs
gene_var = np.array(adata.X.power(2).mean(axis=0) - np.power(adata.X.mean(axis=0), 2)).flatten()
top_gene_indices = np.argsort(gene_var)[-1000:]
genes = adata.var_names[top_gene_indices]

sub_X = adata.X[:, top_gene_indices].toarray()
lib_size = sub_X.sum(axis=1, keepdims=True)
lib_size[lib_size == 0] = 1.0
norm_X = np.log1p((sub_X / lib_size) * 10000.0)
adata_sub = anndata.AnnData(X=norm_X, obs=adata.obs, var=pd.DataFrame(index=genes))

program_correlations = []

for idx, r in p4_p14_matched.iterrows():
    cname = r['cell_name']
    drug = r['drug']
    
    # Cells in plate 4
    c4 = adata_sub.obs[(adata_sub.obs['cell_name'] == cname) & (adata_sub.obs['drug'] == drug) & (adata_sub.obs['plate'] == 'plate4')]
    # Cells in plate 14
    c14 = adata_sub.obs[(adata_sub.obs['cell_name'] == cname) & (adata_sub.obs['drug'] == drug) & (adata_sub.obs['plate'] == 'plate14')]
    
    ast_4 = c4[c4['cell_state_category'] == 'Adaptive_Survival_Tail']
    dom_4 = c4[c4['cell_state_category'] == 'Dominant_Response']
    ast_14 = c14[c14['cell_state_category'] == 'Adaptive_Survival_Tail']
    dom_14 = c14[c14['cell_state_category'] == 'Dominant_Response']
    
    if len(ast_4) >= 2 and len(dom_4) >= 5 and len(ast_14) >= 2 and len(dom_14) >= 5:
        # Delta vector: mean(AST) - mean(Dominant)
        vec4 = np.array(adata_sub[ast_4.index].X.mean(axis=0) - adata_sub[dom_4.index].X.mean(axis=0)).flatten()
        vec14 = np.array(adata_sub[ast_14.index].X.mean(axis=0) - adata_sub[dom_14.index].X.mean(axis=0)).flatten()
        
        rho, p_val = spearmanr(vec4, vec14)
        pr, p_pr = pearsonr(vec4, vec14)
        program_correlations.append({
            'cell_name': cname,
            'drug': drug,
            'n_ast_plate4': len(ast_4),
            'n_ast_plate14': len(ast_14),
            'program_spearman_rho': float(rho),
            'spearman_p': float(p_val),
            'program_pearson_r': float(pr),
            'pearson_p': float(p_pr)
        })

prog_df = pd.DataFrame(program_correlations)
prog_df.to_csv('tables/benchmark/replicate_program_concordance.csv', index=False)

print(f"Evaluated ART program concordance across {len(prog_df)} conditions with tail cells in both plates.")
print(f"  Median Program Spearman rho = {prog_df['program_spearman_rho'].median():+.4f}")
print(f"  Mean Program Spearman rho   = {prog_df['program_spearman_rho'].mean():+.4f}")
print(f"  Conditions with positive concordance: {sum(prog_df['program_spearman_rho'] > 0)} / {len(prog_df)} ({sum(prog_df['program_spearman_rho'] > 0)/len(prog_df)*100:.1f}%)")
print(prog_df.head(10))
