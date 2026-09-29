import os
import sys
import time
import gzip
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix
import anndata as ad

print("Starting SciPlex2 AnnData build...")
t0 = time.time()

# 1. Load gene annotations
gene_path = 'data_raw/external_acute/sciplex2/GSM4150377_sciPlex2_A549_Transcription_Modulators_gene.annotations.txt.gz'
genes = []
with gzip.open(gene_path, 'rt', encoding='utf-8') as f:
    for line in f:
        parts = line.strip().split('\t')
        genes.append(parts)
df_var = pd.DataFrame(genes, columns=['gene_id', 'gene_short_name'])

# Make unique index
counts = {}
idx = []
for name in df_var['gene_short_name']:
    c = counts.get(name, 0)
    if c == 0:
        idx.append(name)
    else:
        idx.append(f"{name}_{c}")
    counts[name] = c + 1

df_var.index = idx
df_var.index.name = 'gene_symbol'
print(f"Loaded {len(df_var)} genes.")

# 2. Load cell annotations and pData
cell_path = 'data_raw/external_acute/sciplex2/GSM4150377_sciPlex2_A549_Transcription_Modulators_cell.annotations.txt.gz'
pdata_path = 'data_raw/external_acute/sciplex2/GSM4150377_sciPlex2_pData.txt.gz'

cells = []
with gzip.open(cell_path, 'rt', encoding='utf-8') as f:
    for line in f:
        cells.append(line.strip().split('\t')[0])

df_obs = pd.read_csv(pdata_path, sep=' ')
df_obs = df_obs.reindex(cells)
df_obs['Cell'] = cells
df_obs.index = cells
df_obs.index.name = 'cell_barcode'

# Parse top_oligo
extracted = df_obs['top_oligo'].str.extract(r'(?P<drug>[^_]+)_(?P<dose>[^_]+)_(?P<well>[A-Z0-9]+)')
for col in extracted.columns:
    df_obs[col] = extracted[col].values

# Dose as numeric
df_obs['dose_numeric'] = pd.to_numeric(df_obs['dose'], errors='coerce')
# Replicate plate
df_obs['plate'] = df_obs['well'].str[0] # 'A' or 'B'
print(f"Loaded {len(df_obs)} cells. Valid top_oligo: {df_obs['top_oligo'].notna().sum()}")

# 3. Read sparse matrix entries
matrix_path = 'data_raw/external_acute/sciplex2/GSM4150377_sciPlex2_A549_Transcription_Modulators_UMI.count.matrix.gz'
print("Reading UMI count matrix into COO sparse...")

row_indices = []
col_indices = []
data_values = []

n_lines = 0
with gzip.open(matrix_path, 'rt', encoding='utf-8') as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p) >= 3:
            row_indices.append(int(p[0]) - 1) # 0-based gene
            col_indices.append(int(p[1]) - 1) # 0-based cell
            data_values.append(float(p[2]))
            n_lines += 1
            if n_lines % 15000000 == 0:
                print(f"  Parsed {n_lines} lines ({n_lines/60135001*100:.1f}%)...")

print(f"Loaded {n_lines} entries. Creating CSR matrix (cells x genes)...")
mat = coo_matrix((data_values, (col_indices, row_indices)), 
                 shape=(len(df_obs), len(df_var)), dtype=np.float32).tocsr()

del row_indices, col_indices, data_values

print("Constructing AnnData...")
adata = ad.AnnData(X=mat, obs=df_obs, var=df_var)

out_path = 'data_processed/sciplex2_adata.h5ad'
os.makedirs('data_processed', exist_ok=True)
adata.write_h5ad(out_path)
t1 = time.time()
print(f"Successfully saved {out_path}: {os.path.getsize(out_path) / (1024*1024):.2f} MB in {t1-t0:.1f}s!")
