import os
import pandas as pd

os.makedirs("manuscript", exist_ok=True)

# 1. FIGURE_RESULT_MAP.csv
fig_map_data = [
    # Figure 1
    {"figure": "Figure 1", "panel": "1A", "result_id": "RES-004", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "cell distribution", "source_data": "code/figures/generate_all_main_figures.py", "script": "code/figures/generate_all_main_figures.py", "manuscript_section": "Result 1"},
    {"figure": "Figure 1", "panel": "1B", "result_id": "RES-004", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "condition", "source_data": "tables/rare_state/g4_condition_flow.csv", "script": "code/figures/generate_all_main_figures.py", "manuscript_section": "Result 1"},
    {"figure": "Figure 1", "panel": "1C", "result_id": "RES-003, RES-008", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "split-well pair", "source_data": "tables/rare_state/dmso_null_calibration_table.csv", "script": "code/figures/generate_all_main_figures.py", "manuscript_section": "Result 1"},
    {"figure": "Figure 1", "panel": "1D", "result_id": "RES-007", "claim_id": "CLM-003", "dataset": "Tahoe-100M", "statistical_unit": "gene", "source_data": "tables/rare_state/ast_differentially_expressed_genes.csv", "script": "code/figures/generate_all_main_figures.py", "manuscript_section": "Result 1, Result 6"},
    # Figure 2
    {"figure": "Figure 2", "panel": "2A", "result_id": "RES-012, RES-016", "claim_id": "CLM-002, CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "condition (N=263)", "source_data": "figures/source_data/fig2_source_data.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2"},
    {"figure": "Figure 2", "panel": "2B", "result_id": "RES-016", "claim_id": "CLM-002, CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "condition (N=263)", "source_data": "figures/source_data/fig2_source_data.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2"},
    {"figure": "Figure 2", "panel": "2C", "result_id": "RES-016", "claim_id": "CLM-002, CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "condition (N=263)", "source_data": "figures/source_data/fig2_source_data.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2"},
    {"figure": "Figure 2", "panel": "2D", "result_id": "RES-016", "claim_id": "CLM-002, CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "condition (N=263)", "source_data": "figures/source_data/fig2_source_data.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2"},
    {"figure": "Figure 2", "panel": "2E", "result_id": "RES-016", "claim_id": "CLM-002, CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "condition (N=263)", "source_data": "figures/source_data/fig2_source_data.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2"},
    {"figure": "Figure 2", "panel": "2F", "result_id": "RES-016", "claim_id": "CLM-002, CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "condition (N=263)", "source_data": "figures/source_data/fig2_source_data.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2"},
    # Figure 3
    {"figure": "Figure 3", "panel": "3A", "result_id": "RES-017", "claim_id": "CLM-006", "dataset": "Tahoe-100M", "statistical_unit": "condition pair (N=79)", "source_data": "figures/source_data/fig3_replicate_source_data.csv", "script": "code/benchmark/run_p2_replicate_reproducibility.py", "manuscript_section": "Result 3"},
    {"figure": "Figure 3", "panel": "3B", "result_id": "RES-017", "claim_id": "CLM-006", "dataset": "Tahoe-100M", "statistical_unit": "condition pair (N=79)", "source_data": "tables/benchmark/replicate_program_concordance.csv", "script": "code/benchmark/run_p2_program_concordance.py", "manuscript_section": "Result 3"},
    {"figure": "Figure 3", "panel": "3C", "result_id": "RES-018", "claim_id": "CLM-007", "dataset": "Tahoe-100M", "statistical_unit": "dose pair (N=133)", "source_data": "figures/source_data/fig3_dose_source_data.csv", "script": "code/benchmark/run_p3_dose_response_architecture.py", "manuscript_section": "Result 4"},
    {"figure": "Figure 3", "panel": "3D", "result_id": "RES-018", "claim_id": "CLM-007", "dataset": "Tahoe-100M", "statistical_unit": "dose pair (N=133)", "source_data": "figures/source_data/fig3_dose_source_data.csv", "script": "code/benchmark/run_p3_dose_response_architecture.py", "manuscript_section": "Result 4"},
    # Figure 4
    {"figure": "Figure 4", "panel": "4A", "result_id": "RES-021", "claim_id": "CLM-008", "dataset": "Tahoe-100M", "statistical_unit": "variance component (N=265)", "source_data": "figures/source_data/fig4_source_data.csv", "script": "code/benchmark/run_p5_variance_partitioning.py", "manuscript_section": "Result 4"},
    {"figure": "Figure 4", "panel": "4B", "result_id": "RES-021", "claim_id": "CLM-008", "dataset": "Tahoe-100M", "statistical_unit": "drug group (N=5)", "source_data": "tables/astr/astr_all_conditions.csv", "script": "code/benchmark/run_p5_variance_partitioning.py", "manuscript_section": "Result 4"},
    # Figure 5
    {"figure": "Figure 5", "panel": "5A", "result_id": "RES-019", "claim_id": "CLM-006", "dataset": "SciPlex2 (GSM4150377)", "statistical_unit": "condition pair (N=28)", "source_data": "figures/source_data/fig5_sciplex_source_data.csv", "script": "code/external_validation/run_sciplex2_p7_benchmark.py", "manuscript_section": "Result 3, Result 5"},
    {"figure": "Figure 5", "panel": "5B", "result_id": "RES-019", "claim_id": "CLM-007", "dataset": "SciPlex2 (GSM4150377)", "statistical_unit": "dose condition (7 doses)", "source_data": "tables/benchmark/sciplex2_dose_architecture.csv", "script": "code/external_validation/run_sciplex2_p7_benchmark.py", "manuscript_section": "Result 4, Result 5"},
    {"figure": "Figure 5", "panel": "5C", "result_id": "RES-020", "claim_id": "CLM-002, CLM-005", "dataset": "Aissa et al. 2021", "statistical_unit": "replicate partition (N=15)", "source_data": "figures/source_data/fig5_aissa_source_data.csv", "script": "code/external_validation/run_aissa_p7_benchmark.py", "manuscript_section": "Result 5"},
    {"figure": "Figure 5", "panel": "5D", "result_id": "RES-020", "claim_id": "CLM-005", "dataset": "Aissa et al. 2021", "statistical_unit": "treatment condition (N=3)", "source_data": "figures/source_data/fig5_aissa_source_data.csv", "script": "code/external_validation/run_aissa_p7_benchmark.py", "manuscript_section": "Result 5"},
    # Figure 6
    {"figure": "Figure 6", "panel": "6A", "result_id": "RES-016", "claim_id": "CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "metric", "source_data": "figures/source_data/fig6_metric_comp_source.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2, Result 6"},
    {"figure": "Figure 6", "panel": "6B", "result_id": "RES-005, RES-009", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "parameter variant", "source_data": "figures/source_data/fig6_def_stab_source.csv", "script": "code/benchmark/run_p11_definition_stability.py", "manuscript_section": "Result 6"},
    {"figure": "Figure 6", "panel": "6C", "result_id": "RES-008", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "negative control model", "source_data": "tables/benchmark/negative_controls_summary.csv", "script": "code/benchmark/run_p12_negative_controls.py", "manuscript_section": "Result 6"},
    # Figure 7
    {"figure": "Figure 7", "panel": "7A", "result_id": "RES-015", "claim_id": "CLM-004, CLM-009", "dataset": "Tahoe-100M vs GSE275199", "statistical_unit": "single cell (HOP62)", "source_data": "figures/source_data/fig7_hop62_source_data.csv", "script": "code/external_dtp/run_g10_blinded_projection.py", "manuscript_section": "Result 7"},
    {"figure": "Figure 7", "panel": "7B", "result_id": "RES-015", "claim_id": "CLM-004, CLM-009", "dataset": "Tahoe-100M vs GSE275199", "statistical_unit": "matched condition pair (N=6)", "source_data": "tables/external_dtp/tahoe_hop62_cells_dtp_projections.csv", "script": "code/external_dtp/check_null3_details.py", "manuscript_section": "Result 7"},
    {"figure": "Figure 7", "panel": "7C", "result_id": "RES-015", "claim_id": "CLM-004, CLM-009", "dataset": "Tahoe-100M vs GSE275199", "statistical_unit": "condition (N=49)", "source_data": "figures/source_data/fig7_pan_dtp_source_data.csv", "script": "code/external_dtp/run_g10_pan_condition_test.py", "manuscript_section": "Result 7"},
]
pd.DataFrame(fig_map_data).to_csv("manuscript/FIGURE_RESULT_MAP.csv", index=False)
print("Saved manuscript/FIGURE_RESULT_MAP.csv")

# 2. TABLE_RESULT_MAP.csv
table_map_data = [
    {"table": "Supplementary Table 1", "result_id": "RES-004", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "condition", "file_path": "tables/rare_state/g4_condition_flow.csv", "script": "code/rare_state/discover_ast_states.py", "manuscript_section": "Result 1, Methods", "description": "Viability gating and cell attrition filtering across conditions"},
    {"table": "Supplementary Table 2", "result_id": "RES-003, RES-008", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "split-well pair", "file_path": "tables/rare_state/dmso_null_calibration_table.csv", "script": "code/rare_state/run_dmso_null_calibration.py", "manuscript_section": "Result 1, Methods", "description": "DMSO technical drift distribution and empirical null floor calibration"},
    {"table": "Supplementary Table 3", "result_id": "RES-007", "claim_id": "CLM-003", "dataset": "Tahoe-100M", "statistical_unit": "gene", "file_path": "tables/rare_state/ast_differentially_expressed_genes.csv", "script": "code/rare_state/discover_ast_states.py", "manuscript_section": "Result 1, Result 6", "description": "Differentially expressed genes in HOP62 ART cells vs dominant response"},
    {"table": "Supplementary Table 4", "result_id": "RES-012, RES-016", "claim_id": "CLM-002, CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "condition", "file_path": "tables/benchmark/comparator_metric_benchmark.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2", "description": "Per-condition metrics for ARTI and C0-C5 comparators (N=263)"},
    {"table": "Supplementary Table 5", "result_id": "RES-016", "claim_id": "CLM-005", "dataset": "Tahoe-100M", "statistical_unit": "metric", "file_path": "tables/benchmark/metric_comparison.csv", "script": "code/benchmark/run_p1_comparator_benchmark.py", "manuscript_section": "Result 2, Result 6", "description": "Multivariate OLS residual variance and cross-validated ML R2 benchmarks"},
    {"table": "Supplementary Table 6", "result_id": "RES-017", "claim_id": "CLM-006", "dataset": "Tahoe-100M", "statistical_unit": "condition pair", "file_path": "tables/benchmark/replicate_reproducibility_pairs.csv", "script": "code/benchmark/run_p2_replicate_reproducibility.py", "manuscript_section": "Result 3", "description": "Cross-plate replicate pair metrics (Plate 4 vs Plate 14, N=79 pairs)"},
    {"table": "Supplementary Table 7", "result_id": "RES-017", "claim_id": "CLM-006", "dataset": "Tahoe-100M", "statistical_unit": "condition pair", "file_path": "tables/replicates/top_quartile_contingency.csv", "script": "code/benchmark/run_p2_replicate_reproducibility.py", "manuscript_section": "Result 3", "description": "Top-quartile categorical contingency table and Cohen's kappa evaluation"},
    {"table": "Supplementary Table 8", "result_id": "RES-017", "claim_id": "CLM-006", "dataset": "Tahoe-100M", "statistical_unit": "condition pair", "file_path": "tables/benchmark/replicate_program_concordance.csv", "script": "code/benchmark/run_p2_program_concordance.py", "manuscript_section": "Result 3", "description": "Transcriptional program rank correlations across independent plates"},
    {"table": "Supplementary Table 9", "result_id": "RES-018", "claim_id": "CLM-007", "dataset": "Tahoe-100M", "statistical_unit": "dose pair", "file_path": "tables/benchmark/dose_response_architecture.csv", "script": "code/benchmark/run_p3_dose_response_architecture.py", "manuscript_section": "Result 4", "description": "Paired two-concentration comparison (0.05 uM vs 5.0 uM, N=133 pairs)"},
    {"table": "Supplementary Table 10", "result_id": "RES-021", "claim_id": "CLM-008", "dataset": "Tahoe-100M", "statistical_unit": "variance component", "file_path": "tables/benchmark/variance_partitioning_report.csv", "script": "code/benchmark/run_p5_variance_partitioning.py", "manuscript_section": "Result 4", "description": "Type II ANOVA variance partitioning across drug, cell line, plate, and residual"},
    {"table": "Supplementary Table 11", "result_id": "RES-019", "claim_id": "CLM-006, CLM-007", "dataset": "SciPlex2", "statistical_unit": "condition pair", "file_path": "tables/benchmark/sciplex2_replicate_pairs.csv", "script": "code/external_validation/run_sciplex2_p7_benchmark.py", "manuscript_section": "Result 3, Result 5", "description": "SciPlex2 independent cross-plate replicate reproducibility and dose series"},
    {"table": "Supplementary Table 12", "result_id": "RES-020", "claim_id": "CLM-002, CLM-005", "dataset": "Aissa et al. 2021", "statistical_unit": "replicate partition", "file_path": "tables/benchmark/aissa_external_benchmark.csv", "script": "code/external_validation/run_aissa_p7_benchmark.py", "manuscript_section": "Result 5", "description": "Aissa et al. PC9 targeted monotherapy vs combination therapy benchmark"},
    {"table": "Supplementary Table 13", "result_id": "RES-005, RES-009", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "parameter variant", "file_path": "tables/benchmark/definition_stability_summary.csv", "script": "code/benchmark/run_p11_definition_stability.py", "manuscript_section": "Result 6", "description": "Definition stability across 10 quantile cutoffs and dispersion methods"},
    {"table": "Supplementary Table 14", "result_id": "RES-008", "claim_id": "CLM-001", "dataset": "Tahoe-100M", "statistical_unit": "negative control model", "file_path": "tables/benchmark/negative_controls_summary.csv", "script": "code/benchmark/run_p12_negative_controls.py", "manuscript_section": "Result 6", "description": "Empirical false-positive rates and ARTI floor across 4 negative controls"},
    {"table": "Supplementary Table 15", "result_id": "RES-015", "claim_id": "CLM-004, CLM-009", "dataset": "Tahoe-100M vs GSE275199", "statistical_unit": "cell", "file_path": "tables/external_dtp/tahoe_hop62_cells_dtp_projections.csv", "script": "code/external_dtp/run_g10_blinded_projection.py", "manuscript_section": "Result 7", "description": "Single-cell DTP projection scores and Null 3 stress-matched pairs for HOP62"},
    {"table": "Supplementary Table 16", "result_id": "RES-015", "claim_id": "CLM-004, CLM-009", "dataset": "Tahoe-100M vs GSE275199", "statistical_unit": "condition", "file_path": "tables/external_dtp/pan_condition_dtp_projections.csv", "script": "code/external_dtp/run_g10_pan_condition_test.py", "manuscript_section": "Result 7", "description": "Pan-condition DTP screening across 49 drug-cell line conditions"},
]
pd.DataFrame(table_map_data).to_csv("manuscript/TABLE_RESULT_MAP.csv", index=False)
print("Saved manuscript/TABLE_RESULT_MAP.csv")

# 3. CLAIM_SENTENCE_MAP.csv
claim_sentence_data = [
    {
        "sentence_id": "SENT-001",
        "section": "Abstract",
        "draft_claim": "Single-cell chemical perturbation distributions exhibit a reproducible minority-tail dimension that is only weakly explained by conventional population-level perturbation metrics.",
        "claim_id": "CLM-002, CLM-005",
        "result_ids": "RES-016",
        "citation_ids": "Srivatsan2020, Zhang2026",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-002",
        "section": "Abstract",
        "draft_claim": "Across 263 conditions in Tahoe-100M, conventional perturbation metrics explained only 15.16% of Acute Response Tail Index (ARTI) variance, leaving 84.84% independent residual variance.",
        "claim_id": "CLM-005",
        "result_ids": "RES-016",
        "citation_ids": "Zhang2026",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-003",
        "section": "Abstract",
        "draft_claim": "In Tahoe-100M, ARTI exhibited moderate cross-plate replicate reproducibility (Pearson r = 0.4173, ICC = 0.3938), while external evaluation in SciPlex2 demonstrated stronger replicate stability (Pearson r = 0.7509, ICC = 0.7155).",
        "claim_id": "CLM-006",
        "result_ids": "RES-017, RES-019",
        "citation_ids": "Srivatsan2020",
        "allowed_status": "QUALIFIED"
    },
    {
        "sentence_id": "SENT-004",
        "section": "Abstract",
        "draft_claim": "Nonlinear machine learning models explained <= 3.09% of ARTI variance in 5-fold cross-validation, demonstrating that tail heterogeneity cannot be reconstructed by non-linear combinations of bulk response metrics.",
        "claim_id": "CLM-005",
        "result_ids": "RES-016",
        "citation_ids": "Peidli2024",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-005",
        "section": "Abstract",
        "draft_claim": "The prespecified 24h acute response tail did not preferentially map to an authentic 8-day drug-tolerant persister benchmark (blinded projection delta rho = -0.0242, p = 0.8354; 47 of 49 conditions null or negative (49 of 49 null after FDR)), establishing a critical negative boundary.",
        "claim_id": "CLM-004, CLM-009",
        "result_ids": "RES-015",
        "citation_ids": "GSE275199",
        "allowed_status": "QUALIFIED"
    },
    {
        "sentence_id": "SENT-006",
        "section": "Result 1",
        "draft_claim": "Vehicle split wells defined an empirical null technical drift floor of 4.57% (P99 = 1.7457), above which genuine acute response tails could be reliably identified.",
        "claim_id": "CLM-001",
        "result_ids": "RES-003, RES-008",
        "citation_ids": "Zhang2026",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-007",
        "section": "Result 1",
        "draft_claim": "Acute response tail cells exhibited balanced cell-cycle phase distributions and no pure G1 arrest across 183 profiled conditions.",
        "claim_id": "CLM-001",
        "result_ids": "RES-006, RES-010",
        "citation_ids": "Zhang2026",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-008",
        "section": "Result 2",
        "draft_claim": "ARTI showed weak correlation with pseudobulk L2 distance (r = -0.0753), energy distance (r = -0.1573), and mean projection magnitude (r = -0.1806; fraction r = +0.0178).",
        "claim_id": "CLM-002, CLM-005",
        "result_ids": "RES-012, RES-016",
        "citation_ids": "Peidli2024",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-009",
        "section": "Result 2",
        "draft_claim": "The residual ARTI component remained reproducible across independent plates (r = +0.4708, p = 1.19e-5), demonstrating that the unexplained variance reflects stable biological signal.",
        "claim_id": "CLM-005, CLM-006",
        "result_ids": "RES-016, RES-017",
        "citation_ids": "Zhang2026",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-010",
        "section": "Result 3",
        "draft_claim": "Cross-plate top-quartile contingency showed 67.09% overall agreement with a positive directional program concordance across 87.3% of conditions.",
        "claim_id": "CLM-006",
        "result_ids": "RES-017",
        "citation_ids": "Zhang2026",
        "allowed_status": "QUALIFIED"
    },
    {
        "sentence_id": "SENT-011",
        "section": "Result 4",
        "draft_claim": "In Tahoe-100M, escalating drug concentration from 0.05 uM to 5.0 uM modulated ARTI (paired t = +5.21, p = 7.18e-7) while remaining decoupled from changes in bulk response magnitude.",
        "claim_id": "CLM-007",
        "result_ids": "RES-018",
        "citation_ids": "Zhang2026",
        "allowed_status": "QUALIFIED"
    },
    {
        "sentence_id": "SENT-012",
        "section": "Result 4",
        "draft_claim": "Type II ANOVA variance partitioning attributed 22.12% of ARTI variance to drug pharmacology (p = 4.91e-17) and 14.04% to cell line background (p = 0.0250).",
        "claim_id": "CLM-008",
        "result_ids": "RES-021",
        "citation_ids": "DepMap2024",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-013",
        "section": "Result 5",
        "draft_claim": "External benchmarking in SciPlex2 confirmed high replicate reproducibility (r = +0.7509) and 63.53% residual non-redundancy beyond comparator metrics.",
        "claim_id": "CLM-005, CLM-006",
        "result_ids": "RES-019",
        "citation_ids": "Srivatsan2020",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-014",
        "section": "Result 5",
        "draft_claim": "In Aissa et al., targeted monotherapy expanded tail dispersion (ARTI = 1.250), whereas dual combination therapy collapsed tail dispersion (ARTI = -0.081) while generating maximal bulk response shift (C0 = 16.20).",
        "claim_id": "CLM-005",
        "result_ids": "RES-020",
        "citation_ids": "Aissa2021",
        "allowed_status": "QUALIFIED"
    },
    {
        "sentence_id": "SENT-015",
        "section": "Result 6",
        "draft_claim": "One-to-one generic stress matching confirmed that ART cells maintain significant excess dispersion (+5.796, p = 0.0230) beyond matched stress controls.",
        "claim_id": "CLM-001",
        "result_ids": "RES-011",
        "citation_ids": "Zhang2026",
        "allowed_status": "ALLOWED"
    },
    {
        "sentence_id": "SENT-016",
        "section": "Result 7",
        "draft_claim": "The generic stress-matched comparison was nominally shifted in the direction opposite to the persister hypothesis (raw p = 0.0408), but did not remain significant after multiple-testing correction (FDR q = 0.1632).",
        "claim_id": "CLM-004, CLM-009",
        "result_ids": "RES-015",
        "citation_ids": "GSE275199",
        "allowed_status": "QUALIFIED"
    },
    {
        "sentence_id": "SENT-017",
        "section": "Result 7",
        "draft_claim": "Acute metallothionein upregulation (MT2A log2FC = +1.87) in 24h tail cells was extinguished in 8-day drug-tolerant persisters (log2FC = -3.41).",
        "claim_id": "CLM-003, CLM-004",
        "result_ids": "RES-007, RES-014",
        "citation_ids": "GSE275199",
        "allowed_status": "ALLOWED"
    }
]
pd.DataFrame(claim_sentence_data).to_csv("manuscript/CLAIM_SENTENCE_MAP.csv", index=False)
print("Saved manuscript/CLAIM_SENTENCE_MAP.csv")
