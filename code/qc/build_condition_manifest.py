import pandas as pd
import ast

samples = pd.read_parquet('data_raw/tahoe/sample_metadata.parquet')
drugs = pd.read_parquet('data_raw/tahoe/drug_metadata.parquet').set_index('drug')

plate_dmso = {}
for plate, grp in samples.groupby('plate'):
    dmso_ids = grp[grp['drug'] == 'DMSO_TF']['sample'].tolist()
    plate_dmso[plate] = ';'.join(dmso_ids)

records = []
for _, row in samples.iterrows():
    s_id = row['sample']
    p = row['plate']
    d = row['drug']
    str_tuple = row['drugname_drugconc']
    conc = '0.0'
    unit = 'uM'
    try:
        val = ast.literal_eval(str_tuple)
        if isinstance(val, list) and len(val) > 0:
            conc = str(val[0][1])
            if len(val[0]) >= 3:
                unit = str(val[0][2])
    except Exception:
        pass
    
    rep = 'rep1'
    if p == 'plate14':
        rep = 'rep2'
    
    moa = ''
    if d in drugs.index:
        broad = str(drugs.loc[d, 'moa-broad'])
        fine = str(drugs.loc[d, 'moa-fine'])
        moa = broad + ' | ' + fine
    
    dmso_ctrl = plate_dmso.get(p, '')
    mito = float(row.get('mean_pcnt_mito', 0.0))
    
    records.append({
        'condition_id': f'{s_id}_{p}_{d}_{conc}{unit}',
        'cell_line': 'Mosaic_pool_50lines',
        'drug': d,
        'dose': f'{conc} {unit}',
        'plate': p,
        'replicate': rep,
        'cell_n_raw': int(row.get('mean_tscp_count', 0)),
        'cell_n_qc': int(row.get('mean_gene_count', 0)),
        'DMSO_control_id': dmso_ctrl,
        'MOA': moa,
        'availability': 'HuggingFace_Parquet',
        'reliability_status': 'PENDING_G2',
        'notes': f'mean_mito={mito:.2f}'
    })

df_out = pd.DataFrame(records)
df_out.to_csv('CONDITION_MANIFEST.csv', index=False)
print('Wrote', len(df_out), 'sample conditions to CONDITION_MANIFEST.csv')
