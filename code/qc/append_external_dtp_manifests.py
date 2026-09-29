import os
import hashlib
import pandas as pd

def get_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            chunk = f.read(8192 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()

data_manifest = pd.read_csv('DATA_MANIFEST.csv')
cache_manifest = pd.read_csv('CACHE_MANIFEST.csv')

new_files = [
    {
        'file_id': 'gse275199_raw_matrix',
        'path': 'data_raw/external_dtp/gse275199/GSE275199_matrix.mtx.gz',
        'source': 'NCBI GEO GSE275199',
        'type': 'external_dtp_raw_counts',
        'format': 'mtx.gz',
        'record_count': 53831939,
        'description': '10x Genomics scRNA-seq raw sparse count matrix (36,601 genes x 8,748 cells)',
        'date_added': '2026-09-28'
    },
    {
        'file_id': 'gse275199_raw_features',
        'path': 'data_raw/external_dtp/gse275199/GSE275199_features.tsv.gz',
        'source': 'NCBI GEO GSE275199',
        'type': 'external_dtp_features',
        'format': 'tsv.gz',
        'record_count': 36601,
        'description': 'Gene annotations and Ensembl IDs for GSE275199 10x matrix',
        'date_added': '2026-09-28'
    },
    {
        'file_id': 'gse275199_raw_barcodes',
        'path': 'data_raw/external_dtp/gse275199/GSE275199_barcodes.tsv.gz',
        'source': 'NCBI GEO GSE275199',
        'type': 'external_dtp_barcodes',
        'format': 'tsv.gz',
        'record_count': 8748,
        'description': '10x cellular barcodes for GSE275199',
        'date_added': '2026-09-28'
    },
    {
        'file_id': 'gse275199_raw_hashtags',
        'path': 'data_raw/external_dtp/gse275199/GSE275199_Samples_hashtags_barcodesequences.txt.gz',
        'source': 'NCBI GEO GSE275199',
        'type': 'external_dtp_hto_metadata',
        'format': 'txt.gz',
        'record_count': 5,
        'description': 'Cell hashing hashtag barcode sequences for H358 and PC9 conditions',
        'date_added': '2026-09-28'
    },
    {
        'file_id': 'gse275199_processed_adata',
        'path': 'data_processed/external_dtp/gse275199_processed_adata.h5ad',
        'source': 'Processed from NCBI GEO GSE275199',
        'type': 'external_dtp_processed_anndata',
        'format': 'h5ad',
        'record_count': 8748,
        'description': 'Demultiplexed, QC-annotated single-cell AnnData for NCI-H358 (Day 0, Day 8 DTP) and PC9 (Day 0, Day 4, Day 11 DTP)',
        'date_added': '2026-09-28'
    }
]

for rec in new_files:
    p = rec['path']
    if os.path.exists(p):
        rec['size_bytes'] = os.path.getsize(p)
        rec['sha256'] = get_sha256(p)
        fid = rec['file_id']
        if fid not in data_manifest['file_id'].values:
            data_manifest = pd.concat([data_manifest, pd.DataFrame([rec])], ignore_index=True)
            print("Added " + fid + " to DATA_MANIFEST.csv")
            
        cid = 'cache_' + fid
        cache_entry = {
            'cache_id': cid,
            'object_name': os.path.basename(p),
            'source_path': p,
            'filters': 'gse275199_dtp_benchmark',
            'code_hash': 'sha256_' + rec['sha256'][:16],
            'creation_date': '2026-09-28',
            'size_bytes': rec['size_bytes'],
            'checksum': rec['sha256'],
            'description': rec['description']
        }
        if cid not in cache_manifest['cache_id'].values:
            cache_manifest = pd.concat([cache_manifest, pd.DataFrame([cache_entry])], ignore_index=True)
            print("Added " + cid + " to CACHE_MANIFEST.csv")

data_manifest.to_csv('DATA_MANIFEST.csv', index=False)
cache_manifest.to_csv('CACHE_MANIFEST.csv', index=False)
print('Updated DATA_MANIFEST.csv and CACHE_MANIFEST.csv successfully.')
