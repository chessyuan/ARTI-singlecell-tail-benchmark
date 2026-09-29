import scanpy as sc
import anndata
import os
import scipy.io

print('Testing 10x mtx loading...')
adata_gse = sc.read_10x_mtx(
    'data_raw/external_dtp/gse275199',
    prefix='GSE275199_',
    var_names='gene_symbols',
    cache=False
)
print('Successfully loaded GSE275199 AnnData:')
print(adata_gse)
print('Shape:', adata_gse.shape)
