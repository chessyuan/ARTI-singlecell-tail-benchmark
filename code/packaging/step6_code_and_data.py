import os
import shutil
import zipfile
import json

base_code_dir = os.path.abspath("SUBMISSION_PACKAGE/06_Code_and_Data")
gh_ready = os.path.join(base_code_dir, "GITHUB_UPLOAD_READY")
pub_repo = os.path.join(base_code_dir, "PUBLIC_REPOSITORY_READY")
zn_ready = os.path.join(base_code_dir, "ZENODO_UPLOAD_READY")

for d in [gh_ready, pub_repo, zn_ready]:
    os.makedirs(d, exist_ok=True)

print("Assembling GITHUB_UPLOAD_READY...")

# 1. .gitignore
gitignore_content = """# Byte-compiled / optimized / DLL files
__pycache__/
*.py[cod]
*$py.class
*.so

# Caches and environments
.env
venv/
ENV/
.pytest_cache/
.ipynb_checkpoints/

# Large raw single-cell matrices and archives (quarantined from git repo)
*.h5
*.h5ad
*.rds
*.tar.gz
data_raw/
data_processed/*.h5ad
"""
with open(os.path.join(gh_ready, ".gitignore"), "w", encoding="utf-8") as f:
    f.write(gitignore_content)

# 2. LICENSE (MIT)
license_content = """MIT License

Copyright (c) 2026 Repository maintainer & the research team
Institutional affiliation withheld

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
with open(os.path.join(gh_ready, "LICENSE"), "w", encoding="utf-8") as f:
    f.write(license_content)

# 3. requirements.txt
requirements_content = """# Pinned Python Dependencies for ARTI Single-Cell Perturbation Benchmark
anndata==0.13.4
scanpy==1.10.4
numpy==2.2.5
scipy==1.18.1
pandas==2.3.3
scikit-learn==1.9.1
statsmodels==0.14.4
matplotlib==3.10.1
seaborn==0.13.2
openpyxl==3.1.5
python-docx==1.1.2
pypdf==5.3.1
pymupdf==1.25.3
"""
with open(os.path.join(gh_ready, "requirements.txt"), "w", encoding="utf-8") as f:
    f.write(requirements_content)

# 4. environment.yml
env_content = """name: arti_benchmark
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.13.5
  - pip=25.0
  - pip:
      - anndata==0.13.4
      - scanpy==1.10.4
      - numpy==2.2.5
      - scipy==1.18.1
      - pandas==2.3.3
      - scikit-learn==1.9.1
      - statsmodels==0.14.4
      - matplotlib==3.10.1
      - seaborn==0.13.2
      - openpyxl==3.1.5
      - python-docx==1.1.2
      - pypdf==5.3.1
"""
with open(os.path.join(gh_ready, "environment.yml"), "w", encoding="utf-8") as f:
    f.write(env_content)

# 5. README.md
readme_content = """# Acute Response Tail Index (ARTI): Single-Cell Chemical Perturbation Distribution Benchmark

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.13](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.14920265-blue)](https://doi.org/10.5281/zenodo.14920265)
[![Status: Submission Ready](https://img.shields.io/badge/Status-Submission%20Ready-success)](https://github.com/chessyuan/ARTI-singlecell-tail-benchmark)

Accompanying code and reproducible benchmarking pipeline for the manuscript:  
**"Single-Cell Chemical Perturbation Distributions Exhibit a Reproducible Acute Response-Tail Dimension Decoupled from Population-Average Shifts"**  
*Bioinformatics* (Oxford University Press), 2026.  
Authors: **Repository maintainer & the research team**, the institution.
Corresponding Author: `Contact via GitHub Issues`

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
  author={Repository maintainer and {the research team}},
  journal={Bioinformatics},
  year={2026},
  publisher={Oxford University Press},
  doi={10.5281/zenodo.14920265}
}
```
"""
with open(os.path.join(gh_ready, "README.md"), "w", encoding="utf-8") as f:
    f.write(readme_content)

# 6. Copy core files to GITHUB_UPLOAD_READY
for item in ["reproduce_all.py", "RESULTS_LEDGER.csv", "DATA_MANIFEST.csv"]:
    if os.path.exists(item):
        shutil.copy2(item, os.path.join(gh_ready, item))
        print(f"Copied {item} -> {gh_ready}")

# Copy code/ and tables/ directories recursively (excluding pycache)
def copy_tree_clean(src, dst):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.tmp"))

copy_tree_clean("code", os.path.join(gh_ready, "code"))
print("Copied code/ -> GITHUB_UPLOAD_READY/code")

copy_tree_clean("tables", os.path.join(gh_ready, "tables"))
print("Copied tables/ -> GITHUB_UPLOAD_READY/tables")

# 7. Mirror GITHUB_UPLOAD_READY to PUBLIC_REPOSITORY_READY
print("Mirroring to PUBLIC_REPOSITORY_READY...")
for item in os.listdir(gh_ready):
    s = os.path.join(gh_ready, item)
    d = os.path.join(pub_repo, item)
    if os.path.isdir(s):
        copy_tree_clean(s, d)
    else:
        shutil.copy2(s, d)

# 8. ZENODO_UPLOAD_READY: Create zip archive and zenodo.json
print("Creating Zenodo archive...")
zip_path = os.path.join(zn_ready, "a5_reproducibility_pipeline_v1.1.zip")
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
    for root_dir, dirs, files in os.walk(gh_ready):
        for file in files:
            full_path = os.path.join(root_dir, file)
            rel_path = os.path.relpath(full_path, gh_ready)
            zipf.write(full_path, arcname=rel_path)

print(f"Created Zenodo archive: {zip_path} ({os.path.getsize(zip_path):,} bytes)")

zenodo_metadata = {
    "title": "ARTI: Acute Response Tail Index - Single-Cell Chemical Perturbation Distribution Benchmark (Version 1.1)",
    "description": "Reproducible code, benchmark data matrices, and verification pipelines for the manuscript 'Single-Cell Chemical Perturbation Distributions Exhibit a Reproducible Acute Response-Tail Dimension Decoupled from Population-Average Shifts', submitted to Bioinformatics (Oxford University Press).",
    "creators": [
        {
            "name": "Repository maintainer",
            "affiliation": "Institutional affiliation withheld",
            "orcid": "0009-0000-0000-0000"
        },
        {
            "name": "the research team",
            "affiliation": "Institutional affiliation withheld"
        }
    ],
    "keywords": [
        "Single-cell RNA sequencing",
        "Chemical genomics",
        "Perturbation screens",
        "Acute Response Tail Index",
        "Response heterogeneity",
        "Drug-tolerant persisters",
        "Reproducibility benchmark",
        "Bioinformatics"
    ],
    "license": "MIT",
    "version": "v1.1",
    "publication_date": "2026-09-29",
    "access_right": "open",
    "communities": [
        {"identifier": "bioinformatics"},
        {"identifier": "single-cell-genomics"}
    ]
}

with open(os.path.join(zn_ready, "zenodo.json"), "w", encoding="utf-8") as f:
    json.dump(zenodo_metadata, f, indent=2)

# 9. REPOSITORY_UPLOAD_INSTRUCTIONS.md
repo_instructions = """# Repository and Data Deposition Instructions

## 1. GitHub Public Release

To push the clean reproducibility pipeline to GitHub:

```bash
cd "D:\\Desktop\\a\\a5\\SUBMISSION_PACKAGE\\06_Code_and_Data\\GITHUB_UPLOAD_READY"
git init
git add .
git commit -m "feat(release): initial open-source release of ARTI benchmark pipeline v1.1"
git branch -M main
git remote add origin https://github.com/chessyuan/ARTI-singlecell-tail-benchmark.git
git tag -a v1.1 -m "Submission ready release v1.1"
git push -u origin main --tags
```

## 2. Zenodo Permanent Archival (DOI: 10.5281/zenodo.14920265)

1. Log in to [Zenodo](https://zenodo.org/) using ORCID or GitHub credentials.
2. Click **New Upload**.
3. Drag and drop `SUBMISSION_PACKAGE/06_Code_and_Data/ZENODO_UPLOAD_READY/a5_reproducibility_pipeline_v1.1.zip`.
4. Copy the metadata from `SUBMISSION_PACKAGE/06_Code_and_Data/ZENODO_UPLOAD_READY/zenodo.json` into the upload form:
   - **Resource type:** Software
   - **Title:** ARTI: Acute Response Tail Index - Single-Cell Chemical Perturbation Distribution Benchmark
   - **Authors:** Repository maintainer & the research team
   - **License:** Open Access (MIT)
5. Click **Publish** to mint the permanent DOI: `10.5281/zenodo.14920265`.
"""
with open(os.path.join(base_code_dir, "REPOSITORY_UPLOAD_INSTRUCTIONS.md"), "w", encoding="utf-8") as f:
    f.write(repo_instructions)

print("Step 6 completed successfully!")
