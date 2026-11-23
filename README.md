# Acute Response Tail Index (ARTI): Single-Cell Chemical Perturbation Distribution Benchmark

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.14920265-blue)](https://doi.org/10.5281/zenodo.14920265)
[![Status: Submission Ready](https://img.shields.io/badge/Status-Submission%20Ready-success)](https://github.com/chessyuan/ARTI-singlecell-tail-benchmark)

Accompanying code and reproducible benchmarking pipeline for the manuscript:  
**"Single-Cell Chemical Perturbation Distributions Exhibit a Reproducible Acute Response-Tail Dimension Decoupled from Population-Average Shifts"**  
*Bioinformatics* (Oxford University Press), 2026.  
Project questions and reproducibility reports: please use GitHub Issues.

---

## Overview

High-throughput single-cell perturbation screens routinely evaluate drug efficacy by measuring population-average expression shifts (such as pseudobulk Euclidean distance or non-parametric Energy distance), leaving extreme minority-tail behaviors largely uncharacterized. This repository provides the complete, deterministic analytical pipeline for:

1. **Acute Response Tail Index (ARTI)**: A standardized metric quantifying excess within-condition dispersion relative to matched vehicle controls in reduced principal component subspace.
2. **Multi-Comparator Decoupling Benchmark**: Systematic benchmarking across 263 conditions in Tahoe-100M demonstrating that ARTI is decoupled from bulk shifts (Pseudobulk L2 $r = -0.0753$, Energy distance $r = -0.1573$), leaving 84.84% of variance unexplained by conventional metrics.
3. **Replicate Reproducibility & Variance Decomposition**: Cross-plate ICC evaluation (Tahoe ICC = 0.3938; SciPlex2 ICC = 0.7155) and Type II ANOVA attributing 22.12% of variance to drug identity, 14.04% to cell context, and 16.35% to technical plate batch.
4. **Longitudinal Persister (DTP) Negative Boundary**: Blinded single-cell projection onto ground-truth 8-day persisters (GSE275199) falsifying the hypothesis that acute 24-hour response tails generically predict long-term drug-tolerant persister programs (Delta rho = -0.0242, $p = 0.8354$; 49/49 conditions null after FDR correction).

---

## Quickstart & Reproduction

### 1. Installation

```bash
git clone https://github.com/chessyuan/ARTI-singlecell-tail-benchmark.git
cd ARTI-singlecell-tail-benchmark
conda env create -f environment.yml
conda activate arti_benchmark
```

### 2. One-Command Master Clean-Room Reproduction

Execute the full benchmark and verification suite:

```bash
python reproduce_all.py
```

### 3. Verify Statistical Integrity & Forensic Linters

```bash
# Verify all 21 milestone numbers match the canonical ledger
python code/manuscript/verify_m1_m4_consistency.py

# Verify zero forbidden causal overclaims or persister misnomers
python code/manuscript/check_forbidden_claims.py

# Independent recomputation of primary statistics
python code/audit/recompute_primary_statistics.py
```

---

## Directory Structure

```
├── code/
│   ├── benchmark/           # Comparator benchmarks, dose-response, OLS residuals
│   ├── rare_state/          # ART identification, QC gating, DMSO null calibration
│   ├── external_validation/ # SciPlex2 and Aissa et al. replication pipelines
│   ├── external_dtp/        # Blinded longitudinal DTP projection suite (GSE275199)
│   ├── manuscript/          # Consistency checkers and claim linters
│   ├── audit/               # Independent statistical recomputation scripts
│   └── figures/             # High-resolution figure rendering routines
├── tables/                  # Processed condition-level benchmark tables & summaries
├── RESULTS_LEDGER.csv       # Canonical ground-truth results ledger (21 milestone numbers)
├── DATA_MANIFEST.csv        # Comprehensive dataset accessions and provenance
├── requirements.txt         # Pinned Python package dependencies
├── environment.yml          # Conda environment definition
├── reproduce_all.py         # Deterministic master execution pipeline
├── LICENSE                  # Open-source MIT License
└── README.md                # Project documentation
```

---

## Citation

```bibtex
@article{Liu2026_ARTI,
  title={Single-Cell Chemical Perturbation Distributions Exhibit a Reproducible Acute Response-Tail Dimension Decoupled from Population-Average Shifts},
  journal={Bioinformatics},
  year={2026},
  publisher={Oxford University Press},
  doi={10.5281/zenodo.14920265}
}
```
