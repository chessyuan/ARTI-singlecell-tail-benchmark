import pandas as pd
import ast

drugs = pd.read_parquet('data_raw/tahoe/drug_metadata.parquet')
samples = pd.read_parquet('data_raw/tahoe/sample_metadata.parquet')

drug_doses = {}
drug_samples = {}
for _, row in samples.iterrows():
    d = row['drug']
    str_tuple = row['drugname_drugconc']
    conc = None
    try:
        val = ast.literal_eval(str_tuple)
        if isinstance(val, list) and len(val) > 0 and len(val[0]) >= 2:
            conc = val[0][1]
    except Exception:
        pass
    
    if d not in drug_doses:
        drug_doses[d] = set()
        drug_samples[d] = 0
    if conc is not None:
        drug_doses[d].add(conc)
    drug_samples[d] += 1

records = []
for _, row in drugs.iterrows():
    d = row['drug']
    doses = sorted(list(drug_doses.get(d, [])))
    dose_str = ';'.join(map(str, doses)) if doses else 'N/A'
    dose_cnt = len(doses)
    cond_cnt = drug_samples.get(d, 0)
    
    moa_broad = row.get('moa-broad', '')
    moa_fine = row.get('moa-fine', '')
    targets = row.get('targets', '')
    
    cytotoxic_keywords = ['dna synthesis', 'topoisomerase', 'tubulin', 'microtubule', 'alkylating', 'antimetabolite', 'chemotherapy', 'rna synthesis', 'intercalat']
    is_cytotoxic = any(k in str(moa_fine).lower() for k in cytotoxic_keywords)
    is_targeted = not is_cytotoxic and pd.notna(targets) and targets != '' and targets is not None and targets != 'None'
    
    notes = str(row.get('gpt-notes-approval', ''))
    if len(notes) > 100:
        notes = notes[:100] + '...'
    
    records.append({
        'drug': d,
        'compound_id': str(row.get('pubchem_cid', '')),
        'PubChem': str(row.get('pubchem_cid', '')),
        'SMILES': str(row.get('canonical_smiles', '')),
        'MOA': f'{moa_broad} | {moa_fine}',
        'targets': str(targets) if pd.notna(targets) else '',
        'approval_status': str(row.get('human-approved', '')),
        'dose_count': dose_cnt,
        'dose_values': dose_str,
        'cell_lines_n': 50,
        'condition_n': cond_cnt,
        'known_targeted_therapy': 'YES' if is_targeted else 'NO',
        'known_cytotoxic': 'YES' if is_cytotoxic else 'NO',
        'notes': notes
    })

df_out = pd.DataFrame(records).sort_values('drug')
df_out.to_csv('DRUG_MANIFEST.csv', index=False)
print('Wrote', len(df_out), 'drugs to DRUG_MANIFEST.csv')
