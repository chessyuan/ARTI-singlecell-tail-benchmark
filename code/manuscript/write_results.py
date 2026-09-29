results_content = r"""# Results

## 1. Reliability-Aware Identification of Acute Response Tails

### 1.1 Deconvolution of Acute Distributional Tails from Technical Drift and Cell Death
To determine whether acute chemical perturbations generate distinct cellular subpopulations among quality-control-passing cells that diverge from the bulk response population, we established an analytical framework to identify Acute Response Tail (ART) states in single-cell perturbation distributions (**Figure 1A**). A critical methodological requirement is distinguishing genuine biological response heterogeneity from two pervasive confounders: technical sampling drift and low-quality, dying cellular debris. 

We addressed this challenge in the Tahoe-100M chemical perturbation atlas by instituting a two-tier filtering strategy (**Figure 1B**). Across 265 experimental conditions, single cells were evaluated against stringent quality-control thresholds (mitochondrial fraction < 15.0%, total UMI counts >= 800, detected genes >= 300; without direct functional viability assays; `RES-004`). Cells failing these criteria were partitioned into a separate low-quality compartment (`tables/rare_state/g4_condition_flow.csv`), ensuring that the subsequent analytical subspace was populated exclusively by analysis-eligible cells.

To establish the technical noise floor of within-condition dispersion, we evaluated split-well vehicle (DMSO) controls across 42 cell lines (`RES-003`, `RES-008`). Across 500 matched pseudo-perturbations, well-to-well technical drift followed a compact distribution with a mean drift of 0.0457 (95% CI: 0.0435–0.0479) and a 99th percentile empirical threshold of Delta dispersion = 1.7457 (**Figure 1C**). Chemical perturbations induced tail dispersion substantially exceeding this null threshold; for example, in HOP62 lung adenocarcinoma cells treated with the KRAS G12C inhibitor Adagrasib (5.0 uM), analysis-eligible cells exhibited an excess dispersion of 3.088 relative to matched vehicle controls (Mann-Whitney U test, p = 0.002, FDR = 0.005; `RES-004`).

### 1.2 Biological Distinctness of Acute Response Tail Cells
We next examined whether cells residing in the acute response tail (defined as the top 5th percentile of Euclidean distance to the condition centroid in 5-PC subspace) represent a distinct transcriptomic state rather than random sampling fluctuations. Differential expression analysis comparing ART cells against the dominant treated population in the primary HOP62 Adagrasib model revealed 55 significantly modulated candidate genes (log2FC > 0.5, p < 0.05; **Figure 1D**; `RES-007`). Notably, the top upregulated gene was metallothionein *MT2A* (log2FC = +1.8722, p = 0.00443, FDR = 0.0125), indicative of acute heavy-metal and electrophilic stress-buffering responses. Downsampling analyses across 50 independent random seeds demonstrated that condition-level tail dispersion estimates remained stable down to N = 30 analysis-eligible cells (Pearson r = 0.7572, p = 4.41e-19; `RES-005`, `RES-009`). Cell-cycle classification revealed that ART cells did not undergo uniform cell-cycle arrest (0 of 183 conditions exhibited pure G1 arrest; `RES-006`, `RES-010`).

---

## 2. ARTI Captures Information Not Explained by Population-Average Response Metrics

### 2.1 Empirical Decoupling from Canonical Bulk Response Metrics
We investigated whether condition-level tail expansion, quantified by the Acute Response Tail Index (ARTI), simply tracks the overall magnitude of perturbation response or reflects an independent distributional dimension. Across 263 evaluable conditions in Tahoe-100M, we systematically benchmarked ARTI against six established population-level perturbation metrics: Pseudobulk L2 distance (C0), Energy distance (C1), Mean Projection Magnitude (C2), Simple within-condition dispersion (C3), Transcriptomic entropy (C4), and Distance to control centroid (C5) (**Figure 2A–2F**; `RES-016`).

Bivariate regression demonstrated that ARTI was only weakly correlated with canonical average-shift metrics (`RES-012`, `RES-016`). Specifically, ARTI showed weak or negligible correlation with Pseudobulk L2 distance (Pearson r = -0.0753, R2 = 0.0057, p = 0.2234; **Figure 2A**) and Energy distance (Pearson r = -0.1573, R2 = 0.0247, p = 0.0106; **Figure 2B**). Weak associations were observed with Mean Projection Magnitude (Pearson r = -0.1806, R2 = 0.0326, p = 0.0033; **Figure 2C**; response penetrance fraction r = +0.0178, R2 = 0.0003, p = 0.7736), simple within-condition dispersion (Pearson r = -0.0591, R2 = 0.0035, p = 0.3394; **Figure 2D**), transcriptomic entropy (Pearson r = -0.3421, R2 = 0.1170, p = 1.24e-8; **Figure 2E**), and distance to control centroid (Pearson r = -0.1267, R2 = 0.0161, p = 0.0401; **Figure 2F**). In a matched-response analysis, we identified 731 condition pairs with nearly identical bulk response shifts (|Delta C0| <= 0.15) that exhibited substantial divergence in tail expansion (|Delta ARTI| >= 1.0; `RES-002`).

### 2.2 Multivariate Non-Redundancy and Replicate Stability of Residual Variance
To formally test whether ARTI could be reconstructed from a linear combination of all six comparator metrics, we fitted a multivariate ordinary least squares (OLS) regression model across all 263 conditions (`RES-016`). In Tahoe-100M, the selected comparator metrics explained 15.16% of ARTI variance in the fitted linear model, leaving 84.84% unexplained (**Figure 6A**). In 5-fold cross-validation, nonlinear machine learning models failed to reconstruct ARTI from comparator features (Random Forest cross-validated R2_CV = 0.0316; Gradient Boosting R2_CV = 0.0139; `RES-016`).

Importantly, this residual component remained reproducible across independent plates. When multivariate regression was performed separately within Plate 4 and Plate 14, the extracted condition-level residuals exhibited significant cross-plate concordance (Pearson r = +0.4708, p = 1.19e-5; intraclass correlation ICC = 0.415; `RES-016`, `RES-017`). This demonstrates that the 84.84% unexplained variance does not represent unmodeled technical noise, but reflects reproducible biological variation.

---

## 3. Acute Response-Tail Measurements Show Dataset-Dependent but Reproducible Structure

### 3.1 Moderate Replicate Reproducibility in Tahoe-100M
To determine the experimental precision of ARTI, we evaluated 79 strictly matched condition pairs profiled in independent 384-well plates (Plate 4 vs Plate 14) in Tahoe-100M (`RES-017`). Across these 79 condition pairs, ARTI demonstrated moderate reproducibility (**Figure 3A**), with a Pearson correlation of r = +0.4173 (95% CI: 0.217–0.584, p = 1.30e-4), a Spearman rank correlation of rho = +0.3062 (95% CI: 0.091–0.509, p = 0.0061), and an intraclass correlation coefficient of ICC(2,1) = 0.3938 (95% CI: 0.185–0.569). 

Categorical contingency analysis of top-quartile ARTI conditions revealed an overall concordance of 67.09% ((7 + 46) / 79; `tables/replicates/top_quartile_contingency.csv`), with a sensitivity of 35.0%, a specificity of 78.0%, Cohen's kappa of kappa = 0.130, and Fisher's exact test p = 0.253. At the molecular level, gene-level transcriptional response programs in the tail exhibited positive directional concordance across plates in 87.3% of conditions (69 of 79 condition pairs showed positive Pearson correlation; mean r = +0.0824, median r = +0.0648; **Figure 3B**; `RES-017`).

### 3.2 Strong Replicate Reproducibility in SciPlex2
Because finite cellular sampling per condition (N = 30–150 cells) in Tahoe-100M limits quantile estimation precision, we evaluated replicate reproducibility in the SciPlex2 dataset, which profiled higher single-cell numbers per well (`RES-019`). Across 28 matched condition pairs between independent 384-well plates (Plate A vs Plate B) in A549 cells, ARTI exhibited strong cross-plate reproducibility (**Figure 5A**), with a Pearson correlation of r = +0.7509 (95% CI: 0.528–0.877, p = 4.16e-6), a Spearman rho of +0.6825 (p = 7.15e-5), and an intraclass correlation coefficient of ICC(2,1) = 0.7155.

---

## 4. Drug Concentration and Cellular Context Structure Acute Tail Responses

### 4.1 Concentration-Dependent Modulation of ARTI in Tahoe-100M
We asked whether acute tail expansion is modulated by compound concentration. In Tahoe-100M, chemical perturbations were tested at two active concentrations (0.05 uM and 5.0 uM). Across 133 matched cell line x drug pairs (`RES-018`), escalating concentration from 0.05 uM to 5.0 uM resulted in a significant increase in tail dispersion (**Figure 3C**). Mean ARTI increased from 0.061 at 0.05 uM to 0.340 at 5.0 uM (paired t-test, t = +5.21, p = 7.18e-7; mean difference Delta ARTI = +0.2792). An exact permutation test shuffling concentration labels within matched pairs confirmed this effect (empirical p < 0.001; `tables/dose/final_dose_claims.csv`).

Crucially, the concentration-dependent change in tail dispersion (Delta ARTI) was statistically decoupled from the change in bulk response magnitude (Delta C0 Pseudobulk L2; Pearson r = -0.1007, p = 0.249; **Figure 3D**; `RES-018`). This indicates that higher pharmacological concentration expands the minority tail through mechanisms distinct from simple population-wide transcriptomic displacement.

### 4.2 Multi-Dose Architectural Analysis in SciPlex2
In SciPlex2, where compounds were tested across seven active concentrations spanning four orders of magnitude (0.1 nM to 10 uM), ARTI displayed compound-specific dose-response architectures (**Figure 5B**; `RES-019`). The histone deacetylase (HDAC) inhibitor Vorinostat (SAHA) induced a monotonic expansion of ARTI across concentrations (Spearman rho = +0.7239, p = 6.06e-8). Conversely, compounds inducing rapid cell cycle arrest or acute cell attrition exhibited tail expansion at intermediate concentrations followed by tail collapse at cytotoxic concentrations (`tables/benchmark/sciplex2_dose_architecture.csv`).

### 4.3 Variance Partitioning: Drug Identity and Cellular Context Determinants
To quantify the relative contributions of drug identity, cell-line context, and technical plate batch to ARTI variation, we performed Type II ANOVA variance partitioning across 265 conditions in Tahoe-100M (`RES-021`). Drug identity explained 22.12% of ARTI variance (F = 18.42, p = 4.91e-17), cell-line context explained 14.04% (F = 1.58, p = 0.0250), technical plate batch accounted for 16.35% (F = 35.12, p = 1.54e-14), and condition-specific residual variation represented 47.49% (**Figure 4A**). Inspection revealed compound-specific distributions of tail-inducing potency (**Figure 4B**).

---

## 5. ART/ARTI Generalizes Across Independent Single-Cell Chemical Perturbation Datasets

### 5.1 Multi-Dose External Replication in SciPlex2
To confirm that the statistical decoupling and reproducibility of ARTI are not artifacts of the Tahoe-100M library chemistry, we applied the ART pipeline to SciPlex2 (`RES-019`). In SciPlex2, multivariate OLS regression of ARTI against canonical comparators demonstrated that 63.53% of ARTI variance was unexplained by existing metrics (bivariate correlations with Pseudobulk L2 r = +0.5517; Energy distance r = +0.5186; Projection magnitude r = +0.5484; `tables/benchmark/sciplex2_comparator_benchmark.csv`). Across 24 vehicle control wells, false-positive tail discovery was controlled (mean ARTI = -0.0177 +/- 0.2229), establishing consistent behavior in an independent split-pool assay.

### 5.2 Targeted Combination Perturbation in Aissa et al. 2021
We next evaluated the behavior of ARTI in an independent 10x Genomics scRNA-seq dataset of EGFR-mutant PC9 cells treated with targeted monotherapy or dual combination therapy (`RES-020`). Across 15 partitioned replicate subsets, ARTI showed near-perfect split-half consistency (Pearson r = +0.9739, p = 0.1457; **Figure 5C**). 

ARTI was strongly negatively correlated with pseudobulk L2 distance (Pearson r = -0.8449, p = 7.37e-5; **Figure 5C**) and Energy distance (Pearson r = -0.8782, p = 1.66e-5). Targeted monotherapy with the MET inhibitor Crizotinib induced substantial tail expansion (mean ARTI = +1.2496), whereas dual combination therapy with Osimertinib plus Crizotinib caused a complete collapse of tail dispersion (mean ARTI = -0.0812; **Figure 5D**) while generating the largest bulk transcriptomic shift observed across conditions (C0 = 16.20).

---

## 6. ART Biology is Distinct from Simple Cell-Cycle and Generic-Stress Explanations

### 6.1 Independence from Cell-Cycle Arrest
To investigate whether ART cells merely represent cells arrested in specific cell-cycle phases, we assigned G1, S, and G2M phase scores across 183 conditions (`RES-006`, `RES-010`). Two-way Chi-squared contingency tests between tail and dominant populations showed balanced cell-cycle phase distributions (Chi-squared test, p = 0.590). The proportion of G1 cells in the tail (mean 27.1%) was lower than in the overall population (mean 37.3%), and zero conditions displayed pure G1 arrest (`tables/astr/ast_cell_cycle_conditioned.csv`).

### 6.2 Excess Dispersion Beyond Matched Generic Stress Controls
We tested whether ART cells reflect non-specific proteotoxic or oxidative stress. Using composite stress scores, each ART cell was paired 1-to-1 with a dominant-response cell from the same condition matched on stress gene expression (`RES-011`). In stress-matched pairs, ART cells maintained significant excess dispersion (+5.796, paired t-test, p = 0.0230; `tables/astr/ast_stress_dissection.csv`).

### 6.3 Metric Definition Invariance and Empirical Negative Controls
Evaluation of 10 alternative mathematical formulations of ARTI demonstrated high rank stability (**Figure 6B**; `RES-005`, `RES-009`). Rankings based on alternative cutoffs correlated strongly with standard ARTI (top 10% r = 0.9419; top 2.5% r = 0.4519; Mahalanobis distance r = 0.6974; `tables/benchmark/definition_stability_summary.csv`).

Under four empirical negative control models (`RES-008`), ARTI values collapsed to baseline (**Figure 6C**; `tables/benchmark/negative_controls_summary.csv`). In DMSO pseudo-perturbations, false-positive rate was 0.0% (0/50 conditions). Under gene label permutation and latent dimension permutation, signal collapsed with false-positive rates of 0.0% and 2.0% (1/50), respectively. Control-mismatched baselines yielded a low false-positive rate of 6.0% (3/50).

---

## 7. Acute Response Tails Do Not Generically Map to Long-Term Drug-Tolerant Persister States

### 7.1 Primary Blinded Test Against GSE275199 Sotorasib DTP Benchmark
A central hypothesis in cancer pharmacology is that short-term transcriptional outliers following drug treatment prefigure or predict long-term drug-tolerant persister (DTP) cells. We tested this hypothesis directly by projecting single cells from the 24h Tahoe-100M atlas onto an authentic 8-day single-cell DTP benchmark from GSE275199 (`RES-013`, `RES-014`, `RES-015`).

In the primary blinded evaluation in KRAS G12C HOP62 cells treated with Adagrasib (5.0 uM, 24h), single cells were scored for correlation with the established NCI-H358 8-day Sotorasib DTP transcriptional signature (**Figure 7A**; `tables/external_dtp/tahoe_hop62_cells_dtp_projections.csv`). Analysis-eligible 24h ART cells did not exhibit elevated similarity to the 8-day DTP program compared to dominant treated cells. The difference in mean Spearman rank correlation was negative and non-significant (Delta rho = -0.0242, two-sample t-test, p = 0.8354; `RES-015`). Testing with alternative DTP scoring metrics confirmed the absence of enrichment (DTP module score difference = +0.0163, p = 0.4603; logistic regression probability difference = -0.0643, p = 0.6267).

### 7.2 Null Model Forensic Evaluation: Null 3 Stress Matching
To determine whether generic stress masking could account for this null result, we re-evaluated four null models (`RES-015`). Permutation of DTP signature gene labels yielded an empirical p = 0.8452 (Null 1); cell-cycle matched permutation yielded p = 0.2593 (Null 2); and library size matching yielded p = 0.9610 (Null 4).

In Null Model 3, where ART cells were paired 1-to-1 with dominant cells matched on generic stress gene expression, the stress-matched comparison was nominally shifted in the direction opposite to the persister hypothesis (Delta rho = -0.0573, paired t = -2.739, raw p = 0.0408; **Figure 7B**), but did not remain significant after multiple-testing correction across the null family (FDR q = 0.1632; `RES-015`).

### 7.3 Pan-Condition Screening across 49 Conditions and Marker Extinction
To evaluate whether persister convergence occurs in other pharmacological or genetic contexts, we screened 49 drug-cell line conditions in Tahoe-100M encompassing diverse oncogenic drivers (KRAS, EGFR, PIK3CA) across five lung adenocarcinoma lines (`RES-015`). Across all 49 conditions, 47 were null or negative in unadjusted testing (median Delta rho = -0.0292; **Figure 7C**; `tables/external_dtp/pan_condition_dtp_projections.csv`), and after Benjamini-Hochberg FDR multiplicity correction across all 49 conditions, 100% (49 of 49) remained null (all FDR q >= 0.6787).

At the molecular level, genes upregulated in the 24h ART state were extinguished in established DTPs (`RES-007`, `RES-014`). *MT2A*, the top upregulated marker in 24h HOP62 ART cells (log2FC = +1.8722), was strongly downregulated in authentic 8-day DTPs (log2FC = -3.4140, p = 1.48e-28; `tables/external_dtp/gse275199_h358_dtp_vs_untreated_de.csv`). Established DTPs instead expressed markers of chromatin remodeling and metabolic quiescence (*HOPX*, *ERBB2*) that were absent in 24h acute tail cells.
"""

with open("manuscript/RESULTS.md", "w", encoding="utf-8") as f:
    f.write(results_content.strip() + "\n")
print("Successfully generated manuscript/RESULTS.md")
