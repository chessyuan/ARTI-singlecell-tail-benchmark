methods_content = r"""# Methods

## 1. Single-Cell Perturbation Data Resources and Provenance

### 1.1 Primary High-Throughput Screen: Tahoe-100M
The primary chemical perturbation atlas was derived from the Tahoe-100M single-cell chemical genomics resource (HuggingFace `tahoebio/Tahoe-100M`). The complete dataset comprises split-pool combinatorial indexing scRNA-seq (sci-RNA-seq3) profiles across human cancer cell lines exposed to diverse pharmacological agents for 24 hours. The analyzed dataset slice includes 10 raw expression shards (`shard_00000`, `shard_00035`, `shard_00039`, `shard_00625`, `shard_00630`, `shard_02594`, `shard_02701`, `shard_02741`, `shard_02748`) containing 282,240 cells total, alongside comprehensive metadata files (`cell_line_metadata.parquet`, `de_sample.parquet`, `drug_metadata.parquet`, `gene_metadata.parquet`, `sample_metadata.parquet`). After strict quality filtering and matched vehicle curation, the finalized primary analysis matrix (Tahoe MVP) encompassed 26,078 QC-passing cells (without direct functional viability assays) across 265 experimental conditions, spanning 42 human cancer cell lines, 7 therapeutic compounds, and 3 independent 384-well experimental plates (Plate 4, Plate 6, Plate 14).

### 1.2 Independent Multi-Dose Acute Validation: SciPlex2
External validation of acute single-cell perturbation reproducibility and multi-dose architecture was conducted using the published SciPlex2 dataset (NCBI GEO `GSM4150377` / `GSE139944`, harmonized in `scPerturb` as `SrivatsanTrapnell2020_sciplex2`). Profiling was performed using sci-RNA-seq2 in A549 lung adenocarcinoma cells treated for 24 hours with four compounds (Nutlin-3a, Dexamethasone, Vorinostat/SAHA, and BMS-345541) across seven active concentrations (0.1, 0.5, 1, 5, 10, 50, and 100 nM or uM nominal) plus matched vehicle controls (dose 0), across two independent 384-well plates (Plate A and Plate B). The dataset contains 23,733 quality-controlled single cells.

### 1.3 Independent Targeted Combination Therapy Benchmark: Aissa et al. 2021
External evaluation of acute tail behavior in targeted monotherapy versus dual combination therapy was performed using single-cell transcriptomics from Aissa et al. 2021 (NCBI GEO `GSE165659`, Zenodo `7041849`, scPerturb `AissaBenevolenskaya2021.h5ad`). This dataset comprises 119,071 cells profiled by 10x Genomics 3' scRNA-seq in EGFR-mutant PC9 non-small cell lung cancer cells subjected to active therapeutic concentrations of the MET inhibitor Crizotinib, the third-generation EGFR inhibitor Osimertinib, or their dual combination, alongside matched vehicle controls.

### 1.4 Longitudinal Ground-Truth Drug-Tolerant Persister Benchmark: GSE275199
To establish a rigorous, authentic ground-truth benchmark for long-term drug-tolerant persister (DTP) programs, we acquired raw 10x Genomics single-cell sequencing matrices from NCBI GEO accession `GSE275199`. This study tracked lung adenocarcinoma models under prolonged targeted therapy: NCI-H358 (KRAS G12C) exposed to Sotorasib for 8 days (yielding established quiescent DTPs, N = 1,291 cells, alongside untreated Day 0 cells, N = 1,861 cells), and PC9 (EGFR exon 19 del) exposed to Osimertinib across Day 0, Day 4 (transient adaptation), and Day 11 (established DTPs, N = 2,222 cells). Raw matrix files (`GSE275199_matrix.mtx.gz`, `GSE275199_features.tsv.gz`, `GSE275199_barcodes.tsv.gz`, and hashtag metadata `GSE275199_Samples_hashtags_barcodesequences.txt.gz`) were demultiplexed and processed into a unified AnnData object containing 8,748 annotated single cells.

---

## 2. Quality Control, Normalization, and Cell Viability Gating

### 2.1 Viability Deconvolution and QC Thresholds
To ensure that transcriptomic outliers represent intact, analysis-eligible cells without functional viability assays rather than late-stage apoptotic fragments or compromised cells, we instituted a two-tier filtering strategy:
1. **Cell Viability Gating**: Individual cellular profiles were classified as QC-excluded if they exhibited elevated mitochondrial transcript abundance (pct_mito >= 15.0%), low unique molecular identifier (UMI) coverage (total_counts < 800), or low detected gene richness (n_genes < 300).
2. **Condition-Level Cell Count Threshold**: Conditions were required to possess a minimum of N_qc >= 25 analysis-eligible cells (or >= 15 cells for vehicle split wells) to be retained for distributional estimation. Conditions with severe cell depletion (< 15 analysis-eligible cells) were flagged as high-attrition and excluded from primary tail estimation. Across all Tahoe conditions, excluded cells were partitioned into a separate dead/dying compartment, ensuring zero contamination of the analysis-eligible analytical space.

### 2.2 Feature Selection and Library Normalization
For the primary perturbation matrix, gene-level raw count vectors were filtered to identify the top 1,000 highly variable genes (HVGs) based on population variance across the atlas. For each single cell, library size was normalized to 10,000 counts per cell, followed by natural logarithmic transformation with a pseudocount of 1:
$$x'_{ij} = \ln\left(1 + \frac{x_{ij}}{\sum_k x_{ik}} \times 10,000\right)$$
where $x_{ij}$ is the raw count for gene $j$ in cell $i$.

---

## 3. Mathematical Definition of the Acute Response Tail (ART)

### 3.1 Principal Component Subspace Decomposition
Dimensionality reduction was performed using Principal Component Analysis (PCA) computed on the log-normalized HVG space across all cells within each experimental plate or batch. The first 5 principal components (PCs) were retained, capturing the dominant axes of biological variance while attenuating technical sequencing shot noise.

### 3.2 Cell Distance and Tail Identification
For an experimental condition $c$ comprising $N_c$ viable single cells with PC coordinates $z_i \in \mathbb{R}^5$ ($i = 1, \dots, N_c$):
1. The condition centroid was established using the component-wise spatial median:
$$\bar{z}_c = \text{median}_{i \in c}(z_i)$$
2. The response distance of each single cell $i$ relative to its condition centroid was calculated as the Euclidean norm:
$$d_i = \|z_i - \bar{z}_c\|_2 = \sqrt{\sum_{k=1}^5 (z_{ik} - \bar{z}_{ck})^2}$$
3. Within condition $c$, the Acute Response Tail (ART) was defined as the subset of cells residing in the top 5th percentile ($q_{0.95}$) of the within-condition distance distribution:
$$\text{ART}_c = \{i \in c \mid d_i \ge \text{Quantile}_{0.95}(d_c)\}$$
The remaining 95% of cells were designated as the Dominant Response population.

---

## 4. Definition and Computation of the Acute Response Tail Index (ARTI)

### 4.1 Condition Dispersion
Condition-level transcriptional dispersion was computed as the arithmetic mean of cellular Euclidean distances:
$$\bar{d}_c = \frac{1}{N_c} \sum_{i \in c} d_i$$

### 4.2 Matched Vehicle Baseline Subtraction
To decouple perturbation-induced dispersion from baseline genetic heterogeneity and cell line-specific clonal diversity, each treated condition was paired with an exact matched vehicle control (DMSO) profiled within the identical cell line and experimental plate:
$$\bar{d}_{\text{DMSO}(c)} = \frac{1}{N_{\text{DMSO}}} \sum_{j \in \text{DMSO}} \|z_j - \bar{z}_{\text{DMSO}}\|_2$$
The excess dispersion induced by chemical perturbation was calculated as:
$$\Delta\bar{d}_c = \bar{d}_c - \bar{d}_{\text{DMSO}(c)}$$

### 4.3 Standardized ARTI Formulation
The Acute Response Tail Index (ARTI) was defined as the baseline-subtracted delta dispersion scaled by a constant factor of 2.0 (corresponding to the empirical median baseline dispersion across Tahoe cell lines):
$$\text{ARTI}_c = \frac{\bar{d}_c - \bar{d}_{\text{DMSO}(c)}}{2.0} = \frac{\Delta\bar{d}_c}{2.0}$$
Under this formulation:
- $\text{ARTI}_c \approx 0.0$ indicates that the drug-treated condition exhibits within-condition dispersion indistinguishable from the vehicle baseline.
- $\text{ARTI}_c > 0.0$ indicates significant tail expansion beyond baseline phenotypic drift.
- $\text{ARTI}_c < 0.0$ denotes transcriptional canalization or phenotypic compression relative to vehicle controls.

Statistical significance of condition-level tail expansion was evaluated using a one-sided Mann-Whitney U test comparing cellular distances $\{d_i\}_{i \in c}$ against matched vehicle distances $\{d_j\}_{j \in \text{DMSO}}$, followed by Benjamini-Hochberg false discovery rate (FDR) adjustment across all conditions.

---

## 5. Technical Calibration and Empirical DMSO Null Floor

To quantify the contribution of well-to-well technical variability, spatial edge effects, and cell sampling drift to ARTI, we evaluated 500 matched pseudo-perturbations constructed from split-half vehicle (DMSO) wells across 42 cell lines. For each pseudo-perturbation pair, one vehicle well was designated as the pseudo-treated condition and compared against its sister vehicle well. Across 500 iterations, the empirical distribution of null delta dispersion yielded a 99th percentile threshold of $\Delta\bar{d}_{\text{null, P99}} = 1.7457$, establishing a strict empirical floor. Across individual split wells, mean technical drift was 0.0457 (95% CI: 0.0435–0.0479). Significant tail conditions were strictly required to satisfy $\text{FDR} < 0.05$ and excess dispersion exceeding the empirical null threshold.

---

## 6. Implementation of Conventional Population-Level Comparators

To evaluate the informational non-redundancy of ARTI, we implemented six canonical population-level perturbation metrics across all evaluated conditions:
1. **C0: Pseudobulk L2 Distance**: The Euclidean distance between the mean log-normalized expression vectors of treated cells ($\mu_c \in \mathbb{R}^{1000}$) and matched vehicle cells ($\mu_{\text{DMSO}} \in \mathbb{R}^{1000}$):
$$\text{C0}_c = \|\mu_c - \mu_{\text{DMSO}}\|_2$$
2. **C1: Energy Distance (scPerturb)**: The statistical Energy Distance between the empirical single-cell distribution of treated cells $P_c$ and matched vehicle cells $P_{\text{DMSO}}$ in 5-PC space, defined as:
$$\text{C1}_c = 2\,\mathbb{E}[\|X - Y\|_2] - \mathbb{E}[\|X - X'\|_2] - \mathbb{E}[\|Y - Y'\|_2]$$
where $X, X' \overset{\text{iid}}{\sim} P_c$ and $Y, Y' \overset{\text{iid}}{\sim} P_{\text{DMSO}}$.
3. **C2: Mean Projection Magnitude**: The Euclidean norm of the centroid displacement vector from vehicle control to treated condition in 5-PC space:
$$\text{C2}_c = \|\bar{z}_c - \bar{z}_{\text{DMSO}}\|_2$$
4. **C3: Simple Within-Condition Dispersion**: The raw mean Euclidean distance of treated cells to the condition centroid, without vehicle baseline subtraction ($\text{C3}_c = \bar{d}_c$).
5. **C4: Transcriptomic Entropy**: The Shannon entropy of the mean expression profile across the top 1,000 HVGs:
$$\text{C4}_c = -\sum_{k=1}^{1000} p_k \ln(p_k), \quad p_k = \frac{\mu_{ck}}{\sum_m \mu_{cm}}$$
6. **C5: Distance to Control Centroid**: The mean Euclidean distance of treated cells directly measured relative to the matched vehicle centroid:
$$\text{C5}_c = \frac{1}{N_c} \sum_{i \in c} \|z_i - \bar{z}_{\text{DMSO}}\|_2$$

*Clarification on Published Metrics*: We note that Song et al. 2025 formulated a specialized Perturbation Score (PS) based on an optimal-transport linear response manifold. Because that specific formulation was not implemented, we refer strictly to our linear displacement metric as Mean Projection Magnitude (C2) and avoid the abbreviation PS to eliminate methodological ambiguity.

---

## 7. Replicate Reproducibility and Concordance Statistics

### 7.1 Statistical Units for Replicate Analysis
To avoid pseudoreplication, all replicate analyses were conducted strictly at the level of experimental conditions (matched cell line x drug x dose combinations across independent plates). Single cells were treated as nested observations within conditions.

### 7.2 Correlation and Intraclass Correlation Analysis
Replicate concordance was evaluated across 79 strictly matched condition pairs between Plate 4 and Plate 14 in Tahoe-100M, and across 28 matched condition pairs between Plate A and Plate B in SciPlex2. Statistical metrics included:
- **Pearson product-moment correlation** ($r$) with two-sided t-test p-values.
- **Spearman rank correlation** ($\rho$).
- **Intraclass Correlation Coefficient** ($\text{ICC}(2,1)$), computed using a two-way mixed-effects model for single measurements absolute agreement:
$$\text{ICC}(2,1) = \frac{\text{MS}_R - \text{MS}_E}{\text{MS}_R + (k - 1)\text{MS}_E + \frac{k}{n}(\text{MS}_C - \text{MS}_E)}$$
where $\text{MS}_R$ is mean square between conditions, $\text{MS}_C$ is mean square between plates ($k=2$), and $\text{MS}_E$ is residual mean square.
- **Bootstrap 95% Confidence Intervals**: Derived from 1,000 non-parametric condition-level resamplings.

### 7.3 Categorical Top-Quartile Concordance
Conditions were partitioned into top-quartile (highest 25% ARTI) versus non-top-quartile within each plate. Concordance was evaluated using a $2 \times 2$ contingency table reporting overall agreement, sensitivity, specificity, Cohen's kappa coefficient ($\kappa$), and two-sided Fisher's exact test.

### 7.4 Transcriptional Program Concordance
For each replicate condition pair, gene-level differential expression vectors (log2 fold change relative to plate-matched vehicle) were calculated across all 1,000 HVGs. Program concordance was defined as the Pearson correlation coefficient between replicate log2FC vectors.

---

## 8. Non-Redundancy and Residual Modeling

### 8.1 Multivariate Ordinary Least Squares (OLS) Regression
To determine whether ARTI can be explained by existing perturbation metrics, we fitted a multivariate linear regression model across 263 conditions in Tahoe-100M:
$$\text{ARTI}_c = \beta_0 + \sum_{m \in \{\text{C0, C1, C2, C3, C4, C5}\}} \beta_m X_{mc} + \epsilon_c$$
The coefficient of determination ($R^2$) was computed, and the proportion of independent residual variance was defined as:
$$\text{Residual Variance (\%)} = (1 - R^2) \times 100\%$$

### 8.2 Cross-Plate Replicate Concordance of the Residual Component
To verify that the residual variance ($\epsilon_c$) represents reproducible biological variation rather than stochastic measurement noise, we extracted residuals from models fitted independently within Plate 4 and Plate 14, and computed the Pearson correlation of residuals across paired conditions ($r_{\text{residual}}$).

### 8.3 Non-Linear Machine Learning Benchmarking
To test whether non-linear interactions among C0-C5 could reconstruct ARTI, we trained Random Forest regressors (100 trees, max depth 5) and Gradient Boosting regressors (learning rate 0.05, max depth 3) using 5-fold cross-validation. Out-of-fold cross-validated $R^2_{\text{CV}}$ was calculated across test folds.

---

## 9. Concentration Modulation and Variance Partitioning

### 9.1 Concentration Escalation Analysis
In Tahoe-100M, chemical perturbations were profiled at two active concentrations (0.05 uM and 5.0 uM). We evaluated 133 matched cell line x drug pairs profiled at both concentrations. Statistical significance of concentration-dependent tail modulation was assessed using a paired two-tailed t-test comparing $\text{ARTI}_{\text{high}}$ versus $\text{ARTI}_{\text{low}}$. Empirical significance was established using an exact label-permutation test (1,000 permutations shuffling concentration labels within pairs).

### 9.2 Type II ANOVA Variance Partitioning
To partition the sources of variation governing ARTI, an ordinary least squares linear model was specified across 265 conditions:
$$\text{ARTI} \sim \text{Drug} + \text{Cell\_Line} + \text{Plate}$$
Type II sum of squares was computed using the `statsmodels` ANOVA framework to quantify the unique variance attributable to compound pharmacology, cell-line genetic context, and technical plate batch.

---

## 10. Sensitivity Analyses and Negative Controls

### 10.1 Cell-Cycle Phase Conditioning
Single-cell cell-cycle phases (G1, S, G2M) were assigned using canonical marker gene sets. Phase distributions between ART cells and dominant response cells were compared using two-way Chi-squared contingency tests. Conditions were evaluated to detect pure G1 arrest (defined as >95% G1 cells in the tail).

### 10.2 Nearest-Neighbor Generic Stress Matching
To decouple compound-specific tail biology from generic proteotoxic and oxidative stress, each ART cell was paired with a dominant-response cell from the same condition matched 1-to-1 on a composite heat-shock and oxidative stress score (based on *HSPA1A*, *HSP90AA1*, *HMOX1*, *SOD1*, *TXN*). Residual dispersion was evaluated in stress-matched pairs using paired t-tests.

### 10.3 Definition Stability
ARTI was recomputed across 10 alternative mathematical formulations, including varied quantile cutoffs (top 1%, 2.5%, 10%), interquartile range (IQR), median absolute deviation (MAD), and Mahalanobis distance in PCA space. Rank stability was quantified by Pearson correlation with the standard formulation.

---

## 11. Longitudinal Drug-Tolerant Persister Benchmark Testing

### 11.1 Blinded Projection on Ground-Truth DTP Signatures
Using the GSE275199 benchmark, ground-truth DTP marker signatures were extracted by differential expression comparing Day 8 Sotorasib DTP cells against Day 0 untreated cells in NCI-H358, and Day 11 Osimertinib DTP cells against Day 0 untreated cells in PC9 (Wilcoxon rank-sum test, FDR < 0.01, log2FC > 1.0). For each single cell in the 24h Tahoe dataset (HOP62 + Adagrasib 5.0 uM, Plate 4), a persister concordance score was computed as the Spearman rank correlation between the cell's normalized gene expression and the DTP fold-change signature across common informative genes. The primary test evaluated the difference in mean persister correlation between 24h ART cells and dominant treated cells ($\Delta\rho = \rho_{\text{ART}} - \rho_{\text{Dom}}$).

### 11.2 Forensic Evaluation of Null Model Specifications
To guard against false-negative conclusions, four null models were tested:
- **Null 1 (Permutation Null)**: Empirical distribution generated by 1,000 random permutations of DTP signature gene labels.
- **Null 2 (Cell-Cycle Matched Null)**: Permutation restricted to cells within the identical cell-cycle phase.
- **Null 3 (Stress-Matched Null)**: ART cells paired 1-to-1 with dominant cells matched on generic stress gene expression.
- **Null 4 (Library Size Null)**: Permutation restricted to cells matched on total UMI counts.

### 11.3 Pan-Condition Screening
The DTP projection test was extended across 49 drug-cell line conditions in Tahoe-100M encompassing diverse oncogenic driver mutations (KRAS, EGFR, PIK3CA) across five non-small cell lung cancer lines.

---

## 12. Statistical Integrity, Software Environments, and Reproducibility

All analyses were executed under Python 3.13.5 (AMD64) on Windows 11 Pro (`10.0.26200`). Core dependencies included `anndata` 0.13.4, `numpy` 2.2.5, `pandas` 2.3.3, `scipy` 1.18.1, `scikit-learn` 1.9.1, `matplotlib` 3.10.1, and `seaborn` 0.13.2. All pseudo-random number generators were initialized with a master seed of 42. All numerical claims, statistical tests, degrees of freedom, and sample sizes were verified against `RESULTS_LEDGER.csv` and `CLAIM_EVIDENCE_MATRIX.csv`. Complete clean-room reproduction is executable via `python reproduce_all.py`.
"""

with open("manuscript/METHODS.md", "w", encoding="utf-8") as f:
    f.write(methods_content.strip() + "\n")
print("Successfully generated manuscript/METHODS.md")
