import pandas as pd
import numpy as np

df = pd.read_csv('RESULTS_LEDGER.csv')

# Drop any previous duplicate records for RES-008..RES-015
df = df.drop_duplicates(subset=['result_id'], keep='last')

# Ensure all 15 records exist and columns are uniform
new_cols = [
    'external_validation_status', 'match_level', 'DTP_duration', 
    'comparison_group', 'null_model', 'stress_matched', 
    'cellcycle_matched', 'primary_metric', 'secondary_metrics'
]

for col in new_cols:
    if col not in df.columns:
        df[col] = 'N/A'

df.loc[df['result_id'].isin(['RES-001', 'RES-002', 'RES-003', 'RES-004', 'RES-005', 'RES-006', 'RES-007']), 'external_validation_status'] = 'INTERNAL_ONLY'

df = df.sort_values(by='result_id').reset_index(drop=True)
df.to_csv('RESULTS_LEDGER.csv', index=False)
print('Cleaned RESULTS_LEDGER.csv. Total unique records:', len(df))
print(df[['result_id', 'hypothesis', 'status', 'metric', 'estimate', 'p_value']])
