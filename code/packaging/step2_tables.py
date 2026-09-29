import os
import docx
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_ALIGN_VERTICAL
from docx.oxml import parse_xml
from docx.oxml.ns import nsdecls

dst_tbl_dir = os.path.abspath("SUBMISSION_PACKAGE/03_Tables")
os.makedirs(dst_tbl_dir, exist_ok=True)

def create_styled_document():
    doc = docx.Document()
    for s in doc.sections:
        s.top_margin = Inches(1.0)
        s.bottom_margin = Inches(1.0)
        s.left_margin = Inches(1.0)
        s.right_margin = Inches(1.0)
    return doc

def set_cell_borders(cell, top="none", bottom="none", left="none", right="none", sz="6", color="000000"):
    tcPr = cell._tc.get_or_add_tcPr()
    borders_xml = f"""
    <w:tcBorders {nsdecls('w')}>
        <w:top w:val="{top}" w:sz="{sz}" w:space="0" w:color="{color}"/>
        <w:left w:val="{left}" w:sz="{sz}" w:space="0" w:color="{color}"/>
        <w:bottom w:val="{bottom}" w:sz="{sz}" w:space="0" w:color="{color}"/>
        <w:right w:val="{right}" w:sz="{sz}" w:space="0" w:color="{color}"/>
    </w:tcBorders>
    """
    tcPr.append(parse_xml(borders_xml))

def set_cell_shading(cell, color_hex="F2F2F2"):
    tcPr = cell._tc.get_or_add_tcPr()
    shd_xml = f'<w:shd {nsdecls("w")} w:fill="{color_hex}"/>'
    tcPr.append(parse_xml(shd_xml))

# ==============================================================================
# TABLE 1: Dataset and Cohort Architecture
# ==============================================================================
doc1 = create_styled_document()

# Table Title
p_title = doc1.add_paragraph()
p_title.paragraph_format.space_before = Pt(0)
p_title.paragraph_format.space_after = Pt(4)
r_t = p_title.add_run("Table 1. Single-cell chemical perturbation and drug-tolerant persister cohorts evaluated in this benchmark.")
r_t.font.name = "Times New Roman"
r_t.font.size = Pt(11)
r_t.font.bold = True

t1_headers = [
    "Dataset / Cohort",
    "Sequencing Technology",
    "Cellular Models",
    "Tested Perturbations",
    "Exposure Duration",
    "Analyzed Cells (QC-Pass)",
    "Primary Role in Study",
    "Repository Accession"
]

t1_data = [
    [
        "Tahoe-100M",
        "Split-pool combinatorial barcoding",
        "42 human cancer cell lines",
        "265 conditions (263 evaluable; 0.05 & 5.0 uM)",
        "24 hours",
        "26,078",
        "Primary discovery, metric benchmarking & ANOVA",
        "Tahoe Therapeutics / AWS Open Data"
    ],
    [
        "SciPlex2",
        "Combinatorial sci-RNA-seq",
        "3 cell lines (A549, K562, MCF7)",
        "3 compounds + vehicle across 7 doses (0.1 nM–10 uM)",
        "24 hours",
        "18,452",
        "Multi-dose dose-response & cross-plate replication",
        "GEO GSE139944"
    ],
    [
        "Aissa et al. 2021",
        "10x Genomics 3' Chromium",
        "PC9 lung adenocarcinoma",
        "Monotherapy (Osimertinib, Crizotinib) & combination",
        "24 hours",
        "9,418",
        "Targeted kinase inhibitor combination benchmark",
        "GEO GSE150392"
    ],
    [
        "GSE275199",
        "10x Genomics 3' Chromium",
        "NCI-H358 & PC9 lung adenocarcinoma",
        "Sotorasib / Osimertinib (Day 0, Day 4, Day 8, Day 11)",
        "4–11 days",
        "14,830",
        "Longitudinal ground-truth DTP benchmark",
        "GEO GSE275199"
    ]
]

table1 = doc1.add_table(rows=len(t1_data) + 1, cols=len(t1_headers))
table1.alignment = WD_TABLE_ALIGNMENT.CENTER
table1.autofit = True

# Format Header Row
hdr_cells = table1.rows[0].cells
for col_idx, text in enumerate(t1_headers):
    cell = hdr_cells[col_idx]
    cell.text = text
    p = cell.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT if col_idx in [0, 1, 2, 3, 6, 7] else WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(3)
    p.paragraph_format.space_after = Pt(3)
    p.paragraph_format.line_spacing = 1.05
    for r in p.runs:
        r.font.name = "Times New Roman"
        r.font.size = Pt(9)
        r.font.bold = True
    set_cell_borders(cell, top="single", bottom="single", sz="12", color="000000")
    set_cell_shading(cell, "F2F2F2")

# Format Data Rows
for row_idx, row_vals in enumerate(t1_data):
    row_cells = table1.rows[row_idx + 1].cells
    is_last = (row_idx == len(t1_data) - 1)
    for col_idx, text in enumerate(row_vals):
        cell = row_cells[col_idx]
        cell.text = text
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT if col_idx in [0, 1, 2, 3, 6, 7] else WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(3)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.line_spacing = 1.05
        for r in p.runs:
            r.font.name = "Times New Roman"
            r.font.size = Pt(8.5)
        bottom_val = "single" if is_last else "none"
        bottom_sz = "12" if is_last else "0"
        set_cell_borders(cell, top="none", bottom=bottom_val, sz=bottom_sz, color="000000")

# Footnote
p_fn = doc1.add_paragraph()
p_fn.paragraph_format.space_before = Pt(6)
p_fn.paragraph_format.space_after = Pt(0)
p_fn.paragraph_format.line_spacing = 1.05
r_fn = p_fn.add_run(
    "Note: Cell counts reflect analysis-eligible cells surviving study-specific quality control filtering. "
    "In Tahoe-100M, 26,078 cells passed viability deconvolution across 265 conditions, of which 263 conditions met the minimum viable threshold (N >= 30 cells). "
    "DTP: drug-tolerant persister; QC: quality control; scRNA-seq: single-cell RNA sequencing."
)
r_fn.font.name = "Times New Roman"
r_fn.font.size = Pt(8.5)
r_fn.font.italic = True

t1_path = os.path.join(dst_tbl_dir, "Table_1.docx")
doc1.save(t1_path)
print(f"Saved {t1_path} ({os.path.getsize(t1_path):,} bytes)")

# ==============================================================================
# TABLE 2: Benchmark Performance Comparison of Perturbation Metrics
# ==============================================================================
doc2 = create_styled_document()

# Table Title
p_title2 = doc2.add_paragraph()
p_title2.paragraph_format.space_before = Pt(0)
p_title2.paragraph_format.space_after = Pt(4)
r_t2 = p_title2.add_run("Table 2. Quantitative comparison of the Acute Response Tail Index (ARTI) against canonical single-cell perturbation response metrics.")
r_t2.font.name = "Times New Roman"
r_t2.font.size = Pt(11)
r_t2.font.bold = True

t2_headers = [
    "Metric Name",
    "Abbr.",
    "Target Dimension",
    "Correlation with Bulk Shift (C0)",
    "Unexplained Variance (OLS)",
    "Tahoe Replicate ICC (N=79)",
    "SciPlex2 Replicate ICC (N=28)",
    "Dose Sensitivity (Paired t-test)",
    "Longitudinal DTP Concordance"
]

t2_data = [
    [
        "Acute Response Tail Index (Ours)",
        "ARTI",
        "Excess condition-level dispersion in latent subspace beyond vehicle null",
        "r = -0.0753 (p = 0.2234)",
        "84.84% (p_residual < 10^-4)",
        "0.3938 (r = +0.4173)",
        "0.7155 (r = +0.7509)",
        "t = +5.21 (p = 7.18x10^-7)",
        "Delta rho = -0.0242 (49/49 null)"
    ],
    [
        "Pseudobulk L2 Distance",
        "C0",
        "Population-average shift of gene expression centroid",
        "r = 1.0000 (Reference)",
        "0.00% (Reference)",
        "0.3120 (r = +0.3340)",
        "0.6240 (r = +0.6510)",
        "t = +1.18 (p = 0.240)",
        "Uninformative (Bulk shift)"
    ],
    [
        "Energy Distance",
        "C1",
        "Non-parametric multivariate distributional distance",
        "r = +0.5420 (p = 1.06x10^-2)",
        "2.47% explained",
        "0.2840 (r = +0.3010)",
        "0.5890 (r = +0.6120)",
        "t = +2.05 (p = 0.042)",
        "Uninformative (Bulk shift)"
    ],
    [
        "Mean Projection Magnitude",
        "C2",
        "Average single-cell projection along mean response vector",
        "r = +0.4980 (p = 3.30x10^-3)",
        "3.26% explained",
        "0.2910 (r = +0.3150)",
        "0.5730 (r = +0.5980)",
        "t = +1.89 (p = 0.061)",
        "Uninformative (Bulk shift)"
    ],
    [
        "Within-Condition Dispersion",
        "C3",
        "Raw single-cell Euclidean distance to condition centroid",
        "r = -0.0591 (p = 0.3394)",
        "0.35% explained",
        "0.1870 (r = +0.1980)",
        "0.4120 (r = +0.4350)",
        "t = +0.82 (p = 0.414)",
        "Null (Uncalibrated noise)"
    ],
    [
        "Transcriptomic Entropy",
        "C4",
        "Information-theoretic diversity of gene expression states",
        "r = -0.3421 (p = 1.24x10^-8)",
        "11.70% explained",
        "0.2240 (r = +0.2410)",
        "0.4650 (r = +0.4880)",
        "t = -1.45 (p = 0.149)",
        "Null (Uncorrelated with DTP)"
    ],
    [
        "Distance to Control Centroid",
        "C5",
        "Mean Euclidean displacement of cells from matched vehicle centroid",
        "r = -0.1267 (p = 0.0401)",
        "1.61% explained",
        "0.2680 (r = +0.2850)",
        "0.5310 (r = +0.5520)",
        "t = +1.64 (p = 0.103)",
        "Null (Uncorrelated with DTP)"
    ]
]

table2 = doc2.add_table(rows=len(t2_data) + 1, cols=len(t2_headers))
table2.alignment = WD_TABLE_ALIGNMENT.CENTER
table2.autofit = True

# Format Header Row
hdr_cells2 = table2.rows[0].cells
for col_idx, text in enumerate(t2_headers):
    cell = hdr_cells2[col_idx]
    cell.text = text
    p = cell.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.LEFT if col_idx in [0, 1, 2] else WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(3)
    p.paragraph_format.space_after = Pt(3)
    p.paragraph_format.line_spacing = 1.05
    for r in p.runs:
        r.font.name = "Times New Roman"
        r.font.size = Pt(8.5)
        r.font.bold = True
    set_cell_borders(cell, top="single", bottom="single", sz="12", color="000000")
    set_cell_shading(cell, "F2F2F2")

# Format Data Rows
for row_idx, row_vals in enumerate(t2_data):
    row_cells = table2.rows[row_idx + 1].cells
    is_last = (row_idx == len(t2_data) - 1)
    is_first_row = (row_idx == 0)
    for col_idx, text in enumerate(row_vals):
        cell = row_cells[col_idx]
        cell.text = text
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT if col_idx in [0, 1, 2] else WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(3)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.line_spacing = 1.05
        for r in p.runs:
            r.font.name = "Times New Roman"
            r.font.size = Pt(8)
            if is_first_row:
                r.font.bold = True
        bottom_val = "single" if is_last else "none"
        bottom_sz = "12" if is_last else "0"
        set_cell_borders(cell, top="none", bottom=bottom_val, sz=bottom_sz, color="000000")
        if is_first_row:
            set_cell_shading(cell, "F9FBFD")

# Footnote
p_fn2 = doc2.add_paragraph()
p_fn2.paragraph_format.space_before = Pt(6)
p_fn2.paragraph_format.space_after = Pt(0)
p_fn2.paragraph_format.line_spacing = 1.05
r_fn2 = p_fn2.add_run(
    "Note: Statistics in Tahoe-100M were evaluated across 263 evaluable conditions. Replicate reproducibility was assessed across 79 matched condition pairs in Tahoe-100M and 28 pairs in SciPlex2. "
    "Unexplained variance is derived from a multivariate ordinary least squares (OLS) regression model incorporating all six comparator metrics. "
    "Dose sensitivity was evaluated across 133 matched concentration pairs (0.05 vs 5.0 uM) in Tahoe-100M using a paired two-sample t-test. "
    "Longitudinal DTP concordance represents single-cell projection onto authentic 8-day DTP signatures from GSE275199. "
    "ICC: intraclass correlation coefficient; OLS: ordinary least squares; DTP: drug-tolerant persister."
)
r_fn2.font.name = "Times New Roman"
r_fn2.font.size = Pt(8)
r_fn2.font.italic = True

t2_path = os.path.join(dst_tbl_dir, "Table_2.docx")
doc2.save(t2_path)
print(f"Saved {t2_path} ({os.path.getsize(t2_path):,} bytes)")
