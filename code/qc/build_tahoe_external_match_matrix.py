import pandas as pd

match_records = [
    {
        'dataset': 'GSE275199',
        'external_cell_line': 'NCI-H358',
        'external_drug': 'Sotorasib (500 nM)',
        'tahoe_matched_cell_line': 'NCI-H358 (CVCL_1559)',
        'tahoe_matched_drug': 'Adagrasib (0.05, 0.5, 5.0 uM)',
        'driver': 'KRAS p.G12C',
        'cancer_type': 'Lung adenocarcinoma (NSCLC)',
        'match_level': 'Level B',
        'single_cell': 'YES (10x GEX)',
        'long_term_duration': '8 days DTP',
        'lineage_tracing': 'NO',
        'validation_feasibility': 'HIGH (Exact cell line H358 + exact driver G12C + exact target class G12Ci)',
        'notes': 'Best targeted NSCLC persister benchmark'
    },
    {
        'dataset': 'GSE275199',
        'external_cell_line': 'PC9',
        'external_drug': 'Osimertinib (150 nM)',
        'tahoe_matched_cell_line': 'LU65 (CVCL_1392, EGFR L858R) / A549 / NCI-H358',
        'tahoe_matched_drug': 'Osimertinib (mesylate) (0.05, 0.5, 5.0 uM)',
        'driver': 'EGFR mutant vs NSCLC',
        'cancer_type': 'Lung adenocarcinoma (NSCLC)',
        'match_level': 'Level C',
        'single_cell': 'YES (10x GEX)',
        'long_term_duration': '4d (early) and 11d (late DTP)',
        'lineage_tracing': 'NO',
        'validation_feasibility': 'HIGH (Exact drug Osimertinib + same cancer NSCLC; PC9 not in Tahoe)',
        'notes': 'Allows cross-temporal acute vs DTP signature comparison'
    },
    {
        'dataset': 'GSE260499',
        'external_cell_line': 'PC9-C2',
        'external_drug': 'Osimertinib (2 uM)',
        'tahoe_matched_cell_line': 'LU65 / NSCLC panel',
        'tahoe_matched_drug': 'Osimertinib (mesylate)',
        'driver': 'EGFR mutant',
        'cancer_type': 'Lung adenocarcinoma (NSCLC)',
        'match_level': 'Level C',
        'single_cell': 'YES (10x GEX, 20 samples)',
        'long_term_duration': '3d, 9d, 14d, 28d, 90d',
        'lineage_tracing': 'YES (Clonal Watermelon barcodes)',
        'validation_feasibility': 'HIGH (Benchmark for longitudinal DTP gene programs & FAK-YAP axis)',
        'notes': 'High-resolution longitudinal single-cell trajectory'
    },
    {
        'dataset': 'GSE196018',
        'external_cell_line': 'A375',
        'external_drug': 'Dabrafenib (1 uM)',
        'tahoe_matched_cell_line': 'RPMI-7951 / C32 / Hs 695T (BRAF V600E melanoma)',
        'tahoe_matched_drug': 'Dabrafenib (0.05, 0.5, 5.0 uM)',
        'driver': 'BRAF p.V600E',
        'cancer_type': 'Melanoma',
        'match_level': 'Level B',
        'single_cell': 'YES (10x GEX)',
        'long_term_duration': 'Weeks (cycling DTPs)',
        'lineage_tracing': 'YES (Clonal barcodes)',
        'validation_feasibility': 'HIGH (Exact drug Dabrafenib + exact driver BRAF V600E + melanoma lineage)',
        'notes': 'Excellent Level B validation for MAPK pathway targeted therapy'
    },
    {
        'dataset': 'GSE193258',
        'external_cell_line': 'PC9, H1975, HCC827, HCC2935',
        'external_drug': 'Osimertinib (500 nM)',
        'tahoe_matched_cell_line': 'NSCLC panel (27 cell lines)',
        'tahoe_matched_drug': 'Osimertinib (mesylate)',
        'driver': 'EGFR mutant',
        'cancer_type': 'Lung adenocarcinoma (NSCLC)',
        'match_level': 'Level C',
        'single_cell': 'NO (Bulk RNA-seq + ATAC)',
        'long_term_duration': '9-14 days DTP',
        'lineage_tracing': 'NO',
        'validation_feasibility': 'MODERATE (Consensus multiomic persister signature)',
        'notes': 'Used to test chromatin accessibility and multiomic convergence'
    },
    {
        'dataset': 'GSE261194',
        'external_cell_line': 'Patient CTCs (BRx-68, etc.)',
        'external_drug': 'Docetaxel (10 nM)',
        'tahoe_matched_cell_line': 'MCF-7 / BT-474 / BT-20 (Breast)',
        'tahoe_matched_drug': 'Docetaxel (0.05, 0.5, 5.0 uM)',
        'driver': 'Various (Luminal / Basal)',
        'cancer_type': 'Breast cancer',
        'match_level': 'Level C',
        'single_cell': 'YES (Smart-seq2)',
        'long_term_duration': '7-10 days DTP',
        'lineage_tracing': 'NO',
        'validation_feasibility': 'MODERATE (Tests whether chemotherapy mitotic poison AST maps to CTC DTP)',
        'notes': 'P27kip1/polyploidy control in persister CTCs'
    }
]

df_match = pd.DataFrame(match_records)
df_match.to_csv('TAHOE_EXTERNAL_MATCH_MATRIX.csv', index=False)
print('Wrote', len(df_match), 'matches to TAHOE_EXTERNAL_MATCH_MATRIX.csv')
