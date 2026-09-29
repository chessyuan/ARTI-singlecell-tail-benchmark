import os

blueprint_content = r"""# Manuscript Blueprint: Acute Single-Cell Drug-Response Tail Heterogeneity

**Document ID**: MANUSCRIPT-BLUEPRINT-001  
**Project**: A5_Acute_Response_Tail_Heterogeneity  
**Phase**: Phase M2 (Manuscript Engineering Blueprint)  
**Status**: ACTIVE BINDING BLUEPRINT  

---

## 1. Title Message & Central Thesis

### 1.1 One-Sentence Title Message
Single-cell chemical perturbation distributions exhibit a reproducible minority-tail dimension that is empirically decoupled from population-average shifts and does not generically converge on established drug-tolerant persister programs.

### 1.2 Central Scientific Question
Is the population-average perturbation response sufficient to characterize single-cell drug responses, or do acute perturbations induce a distinct, reproducible minority-tail dimension of transcriptional variation?

### 1.3 Central Empirical Answer
No, average response is insufficient. Across 265 chemical perturbation conditions in the Tahoe-100M dataset and independent multi-dose platforms (SciPlex2 and Aissa et al.), acute (24h) chemical exposures generate a reproducible minority-tail dimension (Acute Response Tail, ART; quantified by the Acute Response Tail Index, ARTI). In Tahoe-100M, conventional population-level perturbation metrics (pseudobulk L2 distance, energy distance, mean projection magnitude, dispersion, entropy, and control-centroid distance) explain only 15.16% of ARTI variance, leaving 84.84% unexplained residual variance that remains reproducible across independent experimental plates.

### 1.4 Core Negative Boundary
The acute response tail identified here does not preferentially map to established long-term drug-tolerant persister (DTP) states. In a blinded test against an authentic 8-day Sotorasib DTP benchmark (GSE275199), 24h ART cells exhibited no preferential enrichment for the DTP signature compared to dominant treated cells (Spearman correlation difference = -0.0242, p = 0.8354), and 47 of 49 tested conditions across 5 oncology cell lines were null or negative in unadjusted testing (49 of 49 null after Benjamini-Hochberg FDR correction, median delta rho = -0.0292). These findings establish that acute transcriptomic outliers cannot be interpreted as generic pre-existing or early persister states without dedicated longitudinal tracking.

---

## 2. Story Architecture

The manuscript narrative follows a rigorous, unbroken deductive trajectory:

```
Chemical Drug Perturbation
          ↓
Population-Average Shift (Bulk Mean) + Distributional Minority Tail (ART / ARTI)
          ↓
Technical Null Calibration & Viability Deconvolution (Distinguishing Tail from Technical Drift & Apoptosis)
          ↓
Empirical Non-Redundancy (84.84% Residual Variance Beyond Conventional Comparators)
          ↓
Replicate Reproducibility (Moderate in Tahoe, Strong in SciPlex2, High Program Concordance)
          ↓
Pharmacological & Genetic Architecture (Dose Modulation, Factor Partitioning)
          ↓
External Cross-Platform Generalization (SciPlex2 Multi-Dose Series, Aissa et al. Combination Collapse)
          ↓
Definitive Negative Boundary (Falsification of the Generic DTP Hypothesis)
```

---

## 3. Seven Result Modules

### Module 1: Reliability-Aware Identification of Acute Response Tails
- **Scientific Question**: What is the acute response tail (ART), and how can it be rigorously separated from technical noise, apoptosis, cell-cycle arrest, and baseline well drift?
- **Analysis**: Single-cell viable gating; PCA subspace decomposition (5 PCs, top 1,000 HVGs); condition centroid distance distribution; matched DMSO split-well null calibration; differential gene expression analysis.
- **Key Observation**: DMSO vehicle split wells establish an empirical null technical drift floor of 4.57% (P99 = 1.7457). Chemical perturbations induce significant tail expansion exceeding this technical floor (e.g., HOP62 + Adagrasib delta dispersion = 3.088, p = 0.002). Viable ART cells do not exhibit elevated mitochondrial fractions, pure G1 arrest, or cell-cycle skewing.
- **Interpretation**: Acute drug perturbations generate viable transcriptomic tail outliers that represent genuine biological response heterogeneity rather than technical artifacts or dying cells.
- **Figure**: Figure 1 (Panels A, B, C, D).
- **Result IDs**: RES-001, RES-003, RES-004, RES-005, RES-006, RES-008, RES-009, RES-010, RES-011.

### Module 2: ARTI Captures Information Not Explained by Population-Average Metrics
- **Scientific Question**: Does ARTI provide distinct biological information, or is it a trivial proxy for population-average displacement, unmeasured cell death, or baseline entropy?
- **Analysis**: Pairwise bivariate regressions and multivariate ordinary least squares (OLS) regression against six canonical comparator metrics (C0: Pseudobulk L2, C1: Energy distance, C2: Mean projection magnitude, C3: Simple within-condition dispersion, C4: Transcriptomic entropy, C5: Distance to control centroid); 5-fold cross-validated nonlinear machine learning (Random Forest, Gradient Boosting).
- **Key Observation**: Bivariate correlations are weak (pseudobulk L2 r = -0.0753, R2 = 0.0057; energy distance r = -0.1573, R2 = 0.0247; projection magnitude r = -0.1806, R2 = 0.0326; penetrance fraction r = +0.0178, R2 = 0.0003; dispersion r = -0.0591, R2 = 0.0035; entropy r = -0.3421, R2 = 0.1170; dist to control r = -0.1267, R2 = 0.0161). In multivariate regression across 263 conditions, all six comparators explain only 15.16% of ARTI variance (84.84% independent residual variance). Non-linear models explain <= 3.09% of variance in 5-fold cross-validation.
- **Interpretation**: ARTI captures a distributional dimension of chemical perturbation that is statistically decoupled from and non-redundant with conventional mean-shift metrics.
- **Figure**: Figure 2 (Panels A-F) and Figure 6 (Panel A).
- **Result IDs**: RES-012, RES-016.

### Module 3: Acute Response-Tail Measurements Show Dataset-Dependent but Reproducible Structure
- **Scientific Question**: Are acute response-tail measurements reproducible across independent experimental replicates?
- **Analysis**: Cross-plate replicate concordance between independent 384-well plates (Plate 4 vs Plate 14 in Tahoe-100M; Plate A vs Plate B in SciPlex2); intraclass correlation coefficients (ICC[2,1]); top-quartile contingency analysis; gene-level transcriptional program concordance.
- **Key Observation**: In Tahoe-100M, cross-plate ARTI shows moderate reproducibility (Pearson r = +0.4173 [95% CI: 0.217, 0.584], p = 1.30e-4; ICC = 0.3938). Top-quartile agreement is 67.09%, and 87.3% of condition pairs exhibit positive transcriptional program concordance (mean r = +0.0824). The residual ARTI component itself remains reproducible across plates (r = +0.4708, p = 1.19e-5). In SciPlex2, replicate reproducibility is strong (Pearson r = +0.7509, ICC = 0.7155).
- **Interpretation**: Tail dispersion measurements capture stable biological properties with dataset-dependent precision governed by single-cell sampling depth per well.
- **Figure**: Figure 3 (Panels A, B) and Figure 5 (Panel A).
- **Result IDs**: RES-017, RES-019.

### Module 4: Drug Concentration and Cellular Context Structure Acute Tail Responses
- **Scientific Question**: How do compound concentration and cell-line context influence acute tail expansion?
- **Analysis**: Paired comparison of 0.05 uM vs 5.0 uM concentrations across 133 matched conditions in Tahoe-100M; permutation testing of concentration labels; correlation between delta ARTI and delta bulk response; Type II ANOVA variance partitioning across 265 conditions.
- **Key Observation**: ARTI increases significantly across the 100-fold concentration step (mean ARTI increases from 0.061 to 0.340; paired t = +5.21, p = 7.18e-7; permutation p < 0.001), while remaining decoupled from changes in bulk response magnitude (r = -0.1007, p = 0.249). Variance partitioning attributes 22.12% of ARTI variance to drug identity (p = 4.91e-17), 14.04% to cell line background (p = 0.0250), 16.35% to technical plate batch (p = 1.54e-14), and 47.49% to condition-specific residual variation.
- **Interpretation**: Acute tail expansion is a structured pharmacological and genetic phenotype, modulated by drug concentration and cell line permissiveness.
- **Figure**: Figure 3 (Panels C, D) and Figure 4 (Panels A, B).
- **Result IDs**: RES-018, RES-021.

### Module 5: ART/ARTI Generalizes Across Independent Single-Cell Chemical Perturbation Datasets
- **Scientific Question**: Does the ART framework generalize to independent chemical perturbation platforms and experimental designs?
- **Analysis**: Independent evaluation in SciPlex2 (A549, 4 compounds across 7 doses, 2 plates) and Aissa et al. 2021 (PC9, monotherapy vs combination therapy).
- **Key Observation**: In SciPlex2, ARTI exhibits high replicate reproducibility (r = +0.7509, ICC = 0.7155), 63.53% residual variance beyond existing metrics, and monotonic dose expansion under SAHA (Spearman rho = +0.7239, p = 6.06e-8). In Aissa et al., ARTI is strongly negatively correlated with pseudobulk L2 distance (r = -0.8449, p = 7.37e-5); targeted monotherapy (Crizotinib) expands tail dispersion (ARTI = 1.250), whereas dual combination therapy collapses the tail (ARTI = -0.081) despite generating maximum bulk transcriptomic shift (C0 = 16.20).
- **Interpretation**: Tail dynamics represent a general feature of single-cell chemical perturbations, where targeted stress expands distributional dispersion and effective dual-pathway inhibition collapses cellular heterogeneity.
- **Figure**: Figure 5 (Panels A, B, C, D).
- **Result IDs**: RES-019, RES-020.

### Module 6: ART Biology is Distinct from Simple Cell-Cycle and Generic-Stress Explanations
- **Scientific Question**: Is ART simply a manifestation of common cellular arrest or generic stress pathways?
- **Analysis**: Cell-cycle phase conditioning and Chi-squared independence tests; 1-to-1 nearest-neighbor matching on generic stress signatures (heat shock, proteotoxic, oxidative); definition stability testing across 10 alternative cutoffs and dispersion formulations; empirical negative controls (DMSO pseudo-perturbation, gene label permutation, latent dimension permutation).
- **Key Observation**: ART cells show balanced cell-cycle phase distributions (Chi-square p = 0.590); 0 of 183 conditions exhibit pure G1 arrest. Generic stress matching confirms that ART maintains significant excess dispersion beyond matched stress controls (excess dispersion = +5.796, p = 0.0230). ARTI rankings are robust across cutoffs (top 10% r = 0.942, Mahalanobis r = 0.697). All negative controls show false-positive rates <= 0.060.
- **Interpretation**: ART represents a compound- and context-specific transcriptomic tail state distinct from non-specific stress or cell-cycle arrest.
- **Figure**: Figure 1 (Panel D) and Figure 6 (Panels B, C).
- **Result IDs**: RES-006, RES-010, RES-011.

### Module 7: Acute Response Tails Do Not Generically Map to Long-Term Drug-Tolerant Persister States
- **Scientific Question**: Does the acute (24h) response tail state prefigure or predict established (>8-day) drug-tolerant persister (DTP) states?
- **Analysis**: Blinded transcriptomic projection of 24h HOP62 Adagrasib ART cells onto an authentic 8-day Sotorasib DTP single-cell signature (GSE275199, N = 1,291 DTP cells, N = 1,861 untreated); forensic re-evaluation of generic-stress-matched controls (Null 3); pan-condition screening across 49 drug-cell line conditions in Tahoe-100M.
- **Key Observation**: In the primary blinded test, 24h ART cells showed no enrichment for the DTP program relative to dominant treated cells (Spearman rho difference = -0.0242, p = 0.8354). In Null 3 stress-matched pairs, the ART delta rho was nominally negative (-0.0573, raw p = 0.0408, FDR q = 0.1632). Across 49 screened conditions, 47 were null or negative in unadjusted testing (median delta rho = -0.0292), and 49 of 49 remained null after Benjamini-Hochberg FDR correction (all FDR q >= 0.6787). The acute marker MT2A, upregulated in 24h ART (log2FC = +1.87), was extinguished in 8-day DTPs (log2FC = -3.41).
- **Interpretation**: Acute (24h) tail expansion represents immediate transcriptomic diversification that does not generically converge on the established, chromatin-remodeled long-term persister phenotype.
- **Figure**: Figure 7 (Panels A, B, C).
- **Result IDs**: RES-007, RES-013, RES-014, RES-015.

---

## 4. Policy Constraints & Terminology Enforcement

1. **24h Cells in Tahoe**: Termed exclusively **Acute Response Tail (ART)** cells. NEVER termed persisters, DTPs, or precursors.
2. **Persistence Hypothesis**: Formally categorized as a **Negative Boundary / Falsification Result**.
3. **Tahoe Reproducibility**: Transparently reported as **moderate** ($r = 0.4173$, $\text{ICC} = 0.3938$), alongside SciPlex2 ($r = 0.7509$, $\text{ICC} = 0.7155$).
4. **Tahoe Concentration**: Termed **two-concentration modulation** ($0.05\,\mu\text{M}$ vs $5.0\,\mu\text{M}$), never non-monotonic or continuous dose-response.
5. **Aissa Combination Therapy**: Described as **transcriptomic response magnitude** (pseudobulk L2), never extrapolated to unmeasured cell death assays.
"""

with open("manuscript/MANUSCRIPT_BLUEPRINT.md", "w", encoding="utf-8") as f:
    f.write(blueprint_content.strip() + "\n")
print("Successfully generated manuscript/MANUSCRIPT_BLUEPRINT.md")
