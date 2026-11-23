import os
import shutil
import pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
import docx
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import parse_xml
from docx.oxml.ns import nsdecls
import win32com.client

dst_supp_dir = os.path.abspath("SUBMISSION_PACKAGE/04_Supplementary")
os.makedirs(dst_supp_dir, exist_ok=True)

# ==============================================================================
# 1. BUILD Supplementary_Tables.xlsx
# ==============================================================================
print("Creating Supplementary_Tables.xlsx...")
wb = openpyxl.Workbook()
wb.remove(wb.active) # Remove default blank sheet

supp_tables_meta = [
    ("Table S1_QC_Flow", "tables/rare_state/g4_condition_flow.csv", "Table S1. Quality control filtering flow and cell counts across 265 conditions in Tahoe-100M."),
    ("Table S2_DMSO_Null", "tables/rare_state/dmso_null_calibration_table.csv", "Table S2. Empirical vehicle null calibration across 500 matched DMSO split-well pairs."),
    ("Table S3_HOP62_DEGs", "tables/rare_state/ast_differentially_expressed_genes.csv", "Table S3. Differentially expressed genes in HOP62 ART cells under Adagrasib (5.0 uM, 24 h)."),
    ("Table S4_Comparators", "tables/benchmark/comparator_metric_benchmark.csv", "Table S4. Condition-level benchmarking metrics for ARTI and canonical comparators C0-C5 (N=263)."),
    ("Table S5_OLS_ML_Residuals", "tables/benchmark/metric_comparison.csv", "Table S5. Multivariate OLS regression models, residual variance, and machine learning benchmarks."),
    ("Table S6_Replicate_Pairs", "tables/benchmark/replicate_reproducibility_pairs.csv", "Table S6. Replicate reproducibility metrics across 79 matched condition pairs (Plate 4 vs 14)."),
    ("Table S7_TopQuartile", "tables/replicates/top_quartile_contingency.csv", "Table S7. Contingency analysis for top-quartile ARTI conditions across replicate plates."),
    ("Table S8_Program_Concord", "tables/benchmark/replicate_program_concordance.csv", "Table S8. Transcriptional program rank correlations across independent replicate plates."),
    ("Table S9_Dose_Architecture", "tables/benchmark/dose_response_architecture.csv", "Table S9. Paired concentration modulation analysis (0.05 uM vs 5.0 uM, N=133 pairs)."),
    ("Table S10_Dose_Permutation", "tables/dose/final_dose_claims.csv", "Table S10. Exact label-permutation test results for concentration-dependent tail modulation."),
    ("Table S11_SciPlex2_Dose", "tables/benchmark/sciplex2_dose_architecture.csv", "Table S11. Multi-dose dose-response architectures across seven concentrations in SciPlex2."),
    ("Table S12_ANOVA_Variance", "tables/benchmark/variance_partitioning_report.csv", "Table S12. Type II ANOVA variance partitioning across drug, cell line, plate, and residual."),
    ("Table S13_SciPlex2_Comp", "tables/benchmark/sciplex2_comparator_benchmark.csv", "Table S13. Comparator benchmarking, residual variance, and null calibration in SciPlex2."),
    ("Table S14_Aissa_Bench", "tables/benchmark/aissa_external_benchmark.csv", "Table S14. ARTI split-half consistency and combination therapy dynamics in Aissa et al. 2021."),
    ("Table S15_DTP_Projections", "tables/external_dtp/tahoe_hop62_cells_dtp_projections.csv", "Table S15. Blinded single-cell persister projection scores and Null 3 stress-matched pairs."),
    ("Table S16_Pan_DTP_Screen", "tables/external_dtp/pan_condition_dtp_projections.csv", "Table S16. Pan-condition persister projection screening across 49 conditions in Tahoe-100M.")
]

header_fill = PatternFill(start_color="1F497D", end_color="1F497D", fill_type="solid")
header_font = Font(name="Calibri", size=10, bold=True, color="FFFFFF")
title_font = Font(name="Calibri", size=11, bold=True, color="1F497D")
data_font = Font(name="Calibri", size=9.5)
thin_border = Border(
    left=Side(style='thin', color='D9D9D9'),
    right=Side(style='thin', color='D9D9D9'),
    top=Side(style='thin', color='D9D9D9'),
    bottom=Side(style='thin', color='D9D9D9')
)

for tab_name, csv_path, table_title in supp_tables_meta:
    df = pd.read_csv(csv_path)
    ws = wb.create_sheet(title=tab_name[:31]) # Excel sheet name limit 31 chars
    
    # Row 1: Table Title
    ws.cell(row=1, column=1, value=table_title).font = title_font
    
    # Row 3: Headers
    headers = list(df.columns)
    for col_idx, h in enumerate(headers, 1):
        c = ws.cell(row=3, column=col_idx, value=h)
        c.font = header_font
        c.fill = header_fill
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        c.border = thin_border
    
    # Rows 4+: Data
    for row_idx, row_vals in enumerate(df.itertuples(index=False), 4):
        for col_idx, val in enumerate(row_vals, 1):
            c = ws.cell(row=row_idx, column=col_idx)
            if pd.isna(val):
                c.value = "NA"
            elif isinstance(val, (int, float)):
                c.value = val
                if isinstance(val, float):
                    c.number_format = "0.0000" if abs(val) < 100 and abs(val) > 0.0001 else "General"
            else:
                c.value = str(val)
            c.font = data_font
            c.border = thin_border
            c.alignment = Alignment(horizontal="left" if isinstance(val, str) else "right", vertical="center")
            
    # Auto-adjust column widths
    for col in ws.columns:
        col_letter = get_column_letter(col[0].column)
        max_len = 0
        for cell in col[2:]: # ignore row 1 title for width
            if cell.value is not None:
                max_len = max(max_len, len(str(cell.value)))
        ws.column_dimensions[col_letter].width = max(min(max_len + 3, 40), 10)
    
    # Freeze pane at row 4
    ws.freeze_panes = "A4"

xlsx_path = os.path.join(dst_supp_dir, "Supplementary_Tables.xlsx")
wb.save(xlsx_path)
print(f"Saved {xlsx_path} ({os.path.getsize(xlsx_path):,} bytes)")

# ==============================================================================
# 2. BUILD Supplementary_Material.docx & PDF
# ==============================================================================
print("Creating Supplementary_Material.docx...")
doc = docx.Document()
for s in doc.sections:
    s.top_margin = Inches(1.0)
    s.bottom_margin = Inches(1.0)
    s.left_margin = Inches(1.0)
    s.right_margin = Inches(1.0)

# Title
p_t = doc.add_paragraph()
p_t.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p_t.add_run("Supplementary Information")
r.font.name = "Times New Roman"
r.font.size = Pt(18)
r.font.bold = True

p_sub = doc.add_paragraph()
p_sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p_sub.add_run("Single-Cell Chemical Perturbation Distributions Exhibit a Reproducible Acute Response-Tail Dimension Decoupled from Population-Average Shifts\n")
r.font.name = "Times New Roman"
r.font.size = Pt(12)
r.font.bold = True

p_auth = doc.add_paragraph()
p_auth.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p_auth.add_run("Author, affiliation, and contact details omitted from the public repository.\n")
r.font.name = "Times New Roman"
r.font.size = Pt(10.5)

doc.add_paragraph().paragraph_format.space_after = Pt(12)

# Helper function
def add_section_heading(text, level=1):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(16)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.keep_with_next = True
    r = p.add_run(text)
    r.font.name = "Times New Roman"
    r.font.size = Pt(13 if level == 1 else 11.5)
    r.font.bold = True
    return p

def add_body_p(text, bold_prefix="", italic_prefix=""):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.line_spacing = 1.25
    if bold_prefix:
        r_b = p.add_run(bold_prefix)
        r_b.font.name = "Times New Roman"
        r_b.font.size = Pt(10.5)
        r_b.font.bold = True
    if italic_prefix:
        r_i = p.add_run(italic_prefix)
        r_i.font.name = "Times New Roman"
        r_i.font.size = Pt(10.5)
        r_i.font.italic = True
    r = p.add_run(text)
    r.font.name = "Times New Roman"
    r.font.size = Pt(10.5)
    return p

# Supplementary Notes
add_section_heading("Supplementary Note S1: Mathematical Foundations of ARTI and Subspace Geometry", 1)
add_body_p(
    "Conventional high-throughput single-cell perturbation profiling frameworks characterize cellular responses by computing the vector displacement of population means (pseudobulk Euclidean distance, C0) or the non-parametric Wasserstein/Energy distance between multi-cell distributions (C1). While these metrics effectively identify compounds that induce broad transcriptomic shifts across an entire cell population, they are intrinsically insensitive to structured variance expansions occurring in the distributional periphery. Specifically, if a drug shifts the centroid of a cell population by Delta mu while simultaneously expanding the radius of the outer 5th percentile of cells by Delta R, bulk-shift metrics will reflect Delta mu while remaining invariant or underpowered with respect to Delta R. The Acute Response Tail Index (ARTI) was formulated to isolate this peripheral dimension by computing the median Euclidean distance of extreme single cells to the condition's spatial median centroid in an orthogonal 5-dimensional principal component subspace, followed by exact subtraction of the baseline dispersion measured in matched vehicle (DMSO) control wells.",
    bold_prefix="Rationale. "
)
add_body_p(
    "Let X_c in R^(N_c x G) denote the normalized single-cell expression matrix for experimental condition c containing N_c viable cells across G highly variable genes. Principal component projection yields Z_c = X_c V_k in R^(N_c x k), where k = 5. The spatial median centroid mu_c in R^k is determined via the Weiszfeld algorithm to ensure robustness against isolated single-cell artifacts. For each cell i in condition c, the Euclidean distance is d_i = || z_i - mu_c ||_2. The condition dispersion D_c(q) is defined as the median distance among cells falling in the upper tail defined by quantile q (default q = 0.95): D_c(0.95) = median_{i: d_i >= Q_q(d_c)} (d_i). The Acute Response Tail Index is then standardized against matched vehicle controls: ARTI_c = (D_c - D_ctl) / sigma_ctl, where D_ctl and sigma_ctl represent the mean dispersion and standard deviation across vehicle control wells on the identical assay plate.",
    bold_prefix="Mathematical Formulation. "
)

add_section_heading("Supplementary Note S2: Technical Batch Drift Calibration and Statistical Safeguards", 1)
add_body_p(
    "A central vulnerability in single-cell dispersion metrics is the conflation of biological tail expansion with plate-level technical noise, split-well sampling drift, or cell death artifacts. To establish a rigorous statistical floor, we conducted empirical calibration across 500 matched vehicle (DMSO) pseudo-perturbation split-well pairs across 42 distinct cancer cell lines in Tahoe-100M. The empirical null distribution of well-to-well dispersion differences exhibited a mean drift of 0.0457 (95% CI: 0.0435–0.0479) with a 99th-percentile threshold of Delta dispersion = 1.7457. Any condition-level tail expansion failing to exceed this 99th percentile threshold is classified as consistent with technical drift. Furthermore, cell viability was decoupled from tail dispersion by applying strict quality control thresholds: cells with mitochondrial read fractions exceeding 15.0%, UMI counts below 800, or fewer than 300 detected genes were partitioned into a separate low-quality compartment prior to subspace projection.",
    bold_prefix="Empirical Calibration. "
)

add_section_heading("Supplementary Note S3: Type II ANOVA Variance Decomposition Across Biological and Technical Factors", 1)
add_body_p(
    "To quantify the relative contributions of pharmacological mechanism, genetic background, and experimental batch to acute response tail variance, we performed Type II ANOVA across all 265 conditions in the Tahoe-100M atlas. Unlike Type I sequential sums of squares, Type II ANOVA computes the unique variance explained by each factor after conditioning on all other non-interaction terms. The linear model ARTI ~ C(drug) + C(cell_name) + C(plate) attributed 22.12% of total variance to compound identity (F = 18.42, p = 4.91 x 10^-17), 14.04% to cell line context (F = 1.58, p = 0.0250), and 16.35% to technical plate batch (F = 35.12, p = 1.54 x 10^-14), leaving 47.49% residual condition-specific variance. This demonstrates that while technical plate effects exist and must be controlled through matched-control baselines, compound pharmacology exerts the dominant explanatory influence on tail expansion.",
    bold_prefix="Variance Partitioning. "
)

add_section_heading("Supplementary Note S4: Falsification of the Acute-to-Chronic Drug-Tolerant Persister Extrapolation", 1)
add_body_p(
    "A pervasive assumption in cancer single-cell pharmacology is that cells exhibiting extreme transcriptomic divergence within 24 hours of acute drug exposure represent nascent drug-tolerant persisters (DTPs). To empirically test this hypothesis, we performed blinded projection of 24-hour Tahoe-100M single cells onto authentic longitudinal DTP signatures from GSE275199 (encompassing 4-day and 8-day Sotorasib-treated NCI-H358 cells and 11-day Osimertinib-treated PC9 cells). In the primary model (HOP62 under Adagrasib), 24-hour ART cells showed no preferential similarity to the 8-day DTP program compared to the dominant treated population (Delta rho = -0.0242, p = 0.8354). Across a broader screen of 49 conditions encompassing diverse oncogenes (KRAS, EGFR, PIK3CA), 47 of 49 conditions were null or negative in unadjusted testing (median Delta rho = -0.0292), and 49 of 49 conditions were null following Benjamini-Hochberg FDR multiplicity correction. Furthermore, the acute tail marker MT2A (log2FC = +1.8722 in 24-hour ART) was strongly downregulated in Day 8 DTPs (log2FC = -3.4140, p = 1.48 x 10^-28). These findings conclusively demonstrate that acute 24-hour response tails reflect transient, reversible stress-buffering states rather than committed pre-persister lineages.",
    bold_prefix="Longitudinal Falsification. "
)

# Legends for Supplementary Figures
add_section_heading("Supplementary Figure Legends", 1)
add_body_p(
    "Supplementary Figure S1. Empirical DMSO Well-to-Well Technical Drift Distribution. Distribution of pairwise dispersion differences across 500 matched vehicle control split-well pairs across 42 cell lines in Tahoe-100M. The red dashed line denotes the expected null mean (0.00); gray shaded area marks the 99th percentile empirical threshold (Delta dispersion = 1.7457).",
    bold_prefix="Figure S1. "
)
add_body_p(
    "Supplementary Figure S2. Transcriptomic Profiling of HOP62 Acute Response Tail Cells. Volcano plot of differentially expressed genes between acute response tail cells and the dominant response population in HOP62 cells treated with Adagrasib (5.0 uM, 24 h). Red points mark significantly upregulated genes (log2FC > 0.5, p < 0.05), including metallothionein MT2A; blue points mark downregulated genes; gray points indicate non-significant transcripts.",
    bold_prefix="Figure S2. "
)
add_body_p(
    "Supplementary Figure S3. Dissection of Cell-Cycle and Non-Specific Stress Explanations. (A) Cell-cycle phase distribution in acute tail versus dominant populations across 183 conditions, confirming absence of uniform G1 arrest. (B) Nearest-neighbor matching against composite proteotoxic and oxidative stress scores, showing significant excess dispersion beyond matched generic stress controls.",
    bold_prefix="Figure S3. "
)
add_body_p(
    "Supplementary Figure S4. Multi-Comparator Bivariate Scatter Benchmark. Scatter plots displaying bivariate relationships between ARTI and six canonical single-cell perturbation comparators (Pseudobulk L2 C0, Energy distance C1, Mean Projection Magnitude C2, Simple Dispersion C3, Entropy C4, Distance to Control C5) across 263 conditions in Tahoe-100M.",
    bold_prefix="Figure S4. "
)

# Legends for Supplementary Tables
add_section_heading("Supplementary Table Legends", 1)
for tab_name, csv_path, table_title in supp_tables_meta:
    add_body_p(
        f"{table_title} Data matrix available in sheet '{tab_name[:31]}' of Supplementary_Tables.xlsx.",
        bold_prefix=f"{tab_name.split('_')[0]}. "
    )

supp_docx = os.path.join(dst_supp_dir, "Supplementary_Material.docx")
doc.save(supp_docx)
print(f"Saved {supp_docx} ({os.path.getsize(supp_docx):,} bytes)")

# Convert Supplementary_Material.docx to PDF via Word COM
supp_pdf = os.path.join(dst_supp_dir, "Supplementary_Material.pdf")
word = win32com.client.DispatchEx("Word.Application")
word.Visible = False
try:
    doc_com = word.Documents.Open(supp_docx)
    doc_com.SaveAs(supp_pdf, FileFormat=17) # wdFormatPDF
    doc_com.Close()
    print(f"Exported Supplementary_Material.pdf via Word COM: {supp_pdf} ({os.path.getsize(supp_pdf):,} bytes)")
finally:
    word.Quit()

# ==============================================================================
# 3. BUILD Supplementary_Figures.pdf
# ==============================================================================
print("Creating Supplementary_Figures.pdf...")
supp_fig_doc = docx.Document()
for s in supp_fig_doc.sections:
    s.top_margin = Inches(0.8)
    s.bottom_margin = Inches(0.8)
    s.left_margin = Inches(0.8)
    s.right_margin = Inches(0.8)

p_sft = supp_fig_doc.add_paragraph()
p_sft.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = p_sft.add_run("Supplementary Figures S1–S4")
r.font.name = "Times New Roman"
r.font.size = Pt(16)
r.font.bold = True

supp_fig_items = [
    ("Figure S1", "figures/rare_state/fig_dmso_null_calibration.png", "Figure S1. Empirical DMSO Well-to-Well Technical Drift Distribution across 500 matched vehicle split-well pairs in Tahoe-100M."),
    ("Figure S2", "figures/rare_state/fig_ast_volcano_de.png", "Figure S2. Volcano plot of differentially expressed genes in HOP62 ART cells under Adagrasib (5.0 uM, 24 h)."),
    ("Figure S3", "figures/astr/fig_ast_stress_dissection.png", "Figure S3. Dissection of cell-cycle phase conditioning and composite generic stress matching."),
    ("Figure S4", "figures/benchmark/fig_comparator_correlations.png", "Figure S4. Detailed bivariate scatter plots of ARTI against canonical perturbation comparators C0–C5 across 263 conditions.")
]

for label, img_path, caption in supp_fig_items:
    supp_fig_doc.add_page_break()
    p_lbl = supp_fig_doc.add_paragraph()
    p_lbl.paragraph_format.space_before = Pt(0)
    p_lbl.paragraph_format.space_after = Pt(4)
    r_l = p_lbl.add_run(label)
    r_l.font.name = "Times New Roman"
    r_l.font.size = Pt(12)
    r_l.font.bold = True
    
    if os.path.exists(img_path):
        p_img = supp_fig_doc.add_paragraph()
        p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p_img.paragraph_format.space_before = Pt(4)
        p_img.paragraph_format.space_after = Pt(8)
        p_img.add_run().add_picture(os.path.abspath(img_path), width=Inches(6.6))
    
    p_cap = supp_fig_doc.add_paragraph()
    p_cap.paragraph_format.space_before = Pt(4)
    p_cap.paragraph_format.space_after = Pt(12)
    p_cap.paragraph_format.line_spacing = 1.15
    r_c = p_cap.add_run(caption)
    r_c.font.name = "Times New Roman"
    r_c.font.size = Pt(9.5)

supp_fig_docx = os.path.join(dst_supp_dir, "Supplementary_Figures.docx")
supp_fig_pdf = os.path.join(dst_supp_dir, "Supplementary_Figures.pdf")
supp_fig_doc.save(supp_fig_docx)

word = win32com.client.DispatchEx("Word.Application")
word.Visible = False
try:
    doc_com = word.Documents.Open(supp_fig_docx)
    doc_com.SaveAs(supp_fig_pdf, FileFormat=17) # wdFormatPDF
    doc_com.Close()
    print(f"Exported Supplementary_Figures.pdf via Word COM: {supp_fig_pdf} ({os.path.getsize(supp_fig_pdf):,} bytes)")
finally:
    word.Quit()

print("Step 3 completed successfully!")
