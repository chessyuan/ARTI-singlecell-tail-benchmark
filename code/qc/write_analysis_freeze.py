import os

freeze_content = """# A5 ANALYSIS FREEZE RECORD (v1.0-analysis-frozen)

**Document ID**: AUD-FREEZE-001  
**Timestamp**: 2026-09-29 13:05:00 UTC  
**Scope**: Definitive Analysis Freeze for Manuscript Preparation  
**Status**: FROZEN (All Scientific Analysis Locked)  

---

## 1. Git State & Freeze Commit

- **Base Git Commit SHA**: `e226c1e998ca079244bf9dfbd3c1f0535d6aa6c7`
- **Short Commit**: `e226c1e`
- **Commit Message**: `feat(analysis): complete Phase 2/3 benchmarks, audits, tables, and figures prior to analysis freeze`
- **Target Git Tag**: `v1.0-analysis-frozen`
- **Analysis Policy**: DISCOVERY / ANALYSIS MODE officially terminated. MANUSCRIPT ENGINEERING MODE active. No new datasets, no new models, no new hypotheses, and no parameter retuning permitted.

---

## 2. Core Definitions & Implementation Hashes

| Component | Source Script | SHA-256 Hash | Mathematical Definition Summary |
|---|---|---|---|
| **ART Definition** | `code/astr/run_g5_astr_finalization.py` | `f5b25b9b74b366f5cd901ba9b539a9d79f57e7a6b18b63608917a45e0d231a21` | Top 5% Euclidean distance outliers in 5-PC subspace relative to condition centroid |
| **ARTI Implementation** | `code/benchmark/run_p1_comparator_benchmark.py` | `6220edcfcfb95ac2c5f39f19c60bc3bb2e4b3c2d078525df5e52a56ea6dfac9b` | Condition dispersion minus matched DMSO baseline dispersion divided by scale: ARTI = (mean_disp_cond - mean_disp_DMSO) / 2.0 |
| **Reproduction Suite** | `reproduce_all.py` | `b840acb8f82fbb16b14af1af8c1d4cc1a6c7b58be3eda1128cf2687df010c6e7` | 11-step complete clean-room reproduction pipeline |
| **Figure Generator** | `code/figures/generate_all_main_figures.py` | `3ed5af29dec9e13fd1e2de895ea4f6a21dddb205780dae05715226a4a77fcc55` | Generates Figures 1-7 and CSV source data |

---

## 3. Dataset Versions & Provenance

### 3.1 Primary Atlas: Tahoe-100M Slice
- **Source**: HuggingFace `tahoebio/Tahoe-100M`
- **Expression Shards**: 10 parquet shards (`shard_00000`, `shard_00035`, `shard_00039`, `shard_00625`, `shard_00630`, `shard_02594`, `shard_02701`, `shard_02741`, `shard_02748`), 282,240 cells total raw.
- **Metadata Tables**: `cell_line_metadata.parquet`, `de_sample.parquet`, `drug_metadata.parquet`, `gene_metadata.parquet`, `sample_metadata.parquet`.
- **Processed MVP Dataset**: `data_processed/tahoe_mvp_adata.h5ad` (26,078 viable cells across 265 conditions, 42 cell lines, 7 drugs, 3 plates).

### 3.2 External Validation 1: SciPlex2
- **Source**: NCBI GEO `GSM4150377` / `GSE139944`; harmonized in `scPerturb` (`SrivatsanTrapnell2020_sciplex2`).
- **Scope**: 23,733 single cells in A549 across 4 compounds (Nutlin-3a, Dexamethasone, SAHA, BMS-345541) over 7 active doses + vehicle control, profiled across independent 384-well plates (Plate A and Plate B).
- **Processed Artifact**: `data_processed/sciplex2_adata.h5ad`.

### 3.3 External Validation 2: Aissa et al. 2021
- **Source**: Zenodo `7041849` / `scPerturb:AissaBenevolenskaya2021.h5ad` / NCBI GEO `GSE165659`.
- **Scope**: 119,071 cells total; PC9 targeted monotherapies (Crizotinib, Osimertinib) and dual combination therapy vs. vehicle control.
- **Raw File**: `data_raw/external_acute/AissaBenevolenskaya2021.h5ad`.

### 3.4 Longitudinal DTP Ground Truth: GSE275199
- **Source**: NCBI GEO `GSE275199` (10x Genomics scRNA-seq, 36,601 genes x 8,748 cells).
- **Scope**: NCI-H358 Day 0 (N = 1,861) vs. Day 8 Sotorasib DTP (N = 1,291); PC9 Day 0 vs. Day 4 and Day 11 Osimertinib DTP (N = 2,222 Day 11 DTP cells).
- **Raw Files**: `GSE275199_matrix.mtx.gz` (SHA256: `e9c1eda1...`), `GSE275199_features.tsv.gz`, `GSE275199_barcodes.tsv.gz`, `GSE275199_Samples_hashtags_barcodesequences.txt.gz`.
- **Processed Artifact**: `data_processed/external_dtp/gse275199_processed_adata.h5ad` (SHA256: `3bc4f9e9c00ddac761583ef0bbb4700e529404dc21ad0fd39a115d7c22df76cf`).

---

## 4. Manifests & Ledger Hashes

| File | SHA-256 Hash | Records / Scope |
|---|---|---|
| `RESULTS_LEDGER.csv` | `088363ee500b93724573c67990d2726ca161fb25fa160f38560823dbfdd0f9b1` | 21 canonical results (RES-001 through RES-021) |
| `CLAIM_EVIDENCE_MATRIX.csv` | `6fa1dc29878db32de8fdd10a796bd8e28a7bd7b608ca4377c33f0205ce369e8f` | 9 locked claims (CLM-001 through CLM-009) |
| `DATA_MANIFEST.csv` | `830f0e4a14fede87fbe097085f9c0c3064c9711b6c94f8a16b652c56a9311b91` | Complete raw data provenance catalog |
| `EXTERNAL_ACUTE_PERTURBATION_MANIFEST.csv` | `9e38be0dd45ae7d81b78706ac36e1316b93cfa589089410a26c4b99a1a0caedd` | SciPlex2, SciPlex3, PRISM, Aissa metadata |

---

## 5. Main Figure Hashes (figures/main/)

| Figure | File Path | SHA-256 Hash | Description |
|---|---|---|---|
| **Figure 1** | `figures/main/fig1_art_framework_qc_null.png` | `de2f07edf82e4aaf6055a0a11b0952dbbe23665dd52866ff447d2cc93b47c9b0` | Concept, Viability QC flow, DMSO null drift, DEGs |
| **Figure 2** | `figures/main/fig2_arti_vs_average_response.png` | `33a0a7bcc2d0307111ba308c82d8bb50ae4aa9644849e260016d5111154e311c` | Decoupling from C0-C5 comparator metrics (N=263) |
| **Figure 3** | `figures/main/fig3_replicate_dose_reproducibility.png` | `aa85900926b17349cc2b1a1fadb4f6f9bb2008dfc07bf90c4be34ba7b7fc5eba` | Plate 4 vs 14 replicate concordance, dose expansion |
| **Figure 4** | `figures/main/fig4_pan_tahoe_landscape_variance.png` | `d5f4adf374865d1c3fd3f1918364887d791797fa36869d3496d26a341b79dc23` | Type II ANOVA variance partitioning and drug landscape |
| **Figure 5** | `figures/main/fig5_external_acute_validation.png` | `f5953d355604bafd05359b25a191bd83a110980c330795a8e0bcc29bf4da2fde` | SciPlex2 (replicates & dose) + Aissa et al. (combo) |
| **Figure 6** | `figures/main/fig6_metric_benchmark_competition.png` | `3ee0891c129fd156d67f926ee9fdba51bb4d27d20158ccd0aa13118c03006eda` | Residual variance competition, definition stability, negative controls |
| **Figure 7** | `figures/main/fig7_dtp_negative_boundary.png` | `215898f753736c57981d3ae7820fb46b9474a07c0130d3d1925f69c8b861ed10` | Primary blinded DTP test, Null 3 forensic recheck, 49-condition screen |

---

## 6. Computational Environment Specification

- **OS / Platform**: Windows 11 Pro 64-bit (`Windows-11-10.0.26200-SP0`)
- **Python Runtime**: `Python 3.13.5` (64-bit AMD64)
- **Key Scientific Packages**:
  - `anndata`: `0.13.4`
  - `numpy`: `2.2.5`
  - `pandas`: `2.3.3`
  - `scipy`: `1.18.1`
  - `scikit-learn`: `1.9.1`
  - `matplotlib`: `3.10.1`
  - `seaborn`: `0.13.2`

---

## 7. Post-Freeze Amendment Policy

Modifications to any frozen analysis, table, or figure are prohibited unless one of the following five criteria is formally satisfied:
1. **Definite Code Bug**: A demonstrable algorithmic or indexing error discovered in the pipeline.
2. **Source Data Inconsistency**: A numerical discrepancy between the text/table and verified source data.
3. **Statistical Unit Violation**: Identification of pseudoreplication or incorrect degrees of freedom.
4. **Clean-Room Failure**: Failure of `reproduce_all.py` to run to completion from source artifacts.
5. **Critical Literature or Methodological Error**: Falsification of a citation or factual error in methodology.

Any permitted change MUST be formally documented in `audit/POST_FREEZE_AMENDMENT.md` prior to code or text alteration.
"""

os.makedirs('audit', exist_ok=True)
with open('audit/ANALYSIS_FREEZE.md', 'w', encoding='utf-8') as f:
    f.write(freeze_content.strip() + '\n')
print('Successfully generated audit/ANALYSIS_FREEZE.md')
