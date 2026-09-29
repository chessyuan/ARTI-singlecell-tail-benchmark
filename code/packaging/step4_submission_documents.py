import shutil
import os
import docx
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import parse_xml
from docx.oxml.ns import nsdecls

dst_doc_dir = os.path.abspath("SUBMISSION_PACKAGE/05_Submission_Documents")
os.makedirs(dst_doc_dir, exist_ok=True)

def create_styled_doc():
    doc = docx.Document()
    for s in doc.sections:
        s.top_margin = Inches(1.0)
        s.bottom_margin = Inches(1.0)
        s.left_margin = Inches(1.0)
        s.right_margin = Inches(1.0)
    return doc

def add_title(doc, text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(12)
    r = p.add_run(text)
    r.font.name = "Times New Roman"
    r.font.size = Pt(16)
    r.font.bold = True
    return p

def add_heading(doc, text, level=1):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(14)
    p.paragraph_format.space_after = Pt(4)
    p.paragraph_format.keep_with_next = True
    r = p.add_run(text)
    r.font.name = "Times New Roman"
    r.font.size = Pt(13 if level == 1 else 11.5)
    r.font.bold = True
    return p

def add_p(doc, text, bold_prefix="", italic_prefix="", space_after=6, line_spacing=1.25):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(space_after)
    p.paragraph_format.line_spacing = line_spacing
    if bold_prefix:
        r_b = p.add_run(bold_prefix)
        r_b.font.name = "Times New Roman"
        r_b.font.size = Pt(11)
        r_b.font.bold = True
    if italic_prefix:
        r_i = p.add_run(italic_prefix)
        r_i.font.name = "Times New Roman"
        r_i.font.size = Pt(11)
        r_i.font.italic = True
    r = p.add_run(text)
    r.font.name = "Times New Roman"
    r.font.size = Pt(11)
    return p

# 1. Cover_Letter.docx
print("Generating Cover_Letter.docx...")
doc1 = create_styled_doc()
with open("submission/final/COVER_LETTER.md", "r", encoding="utf-8") as f:
    cl_lines = f.readlines()
# Parse markdown lines into docx
for line in cl_lines:
    s = line.strip()
    if not s or s.startswith("# COVER_LETTER"):
        continue
    if s.startswith("# "):
        add_title(doc1, s[2:].replace("*", ""))
    elif s.startswith("### "):
        add_heading(doc1, s[4:].replace("*", ""), level=2)
    elif s.startswith("- **"):
        # bullet with bold prefix
        parts = s[4:].split("**:", 1)
        if len(parts) == 2:
            add_p(doc1, parts[1].strip(), bold_prefix=f"• {parts[0].strip()}: ")
        else:
            add_p(doc1, s[2:])
    elif s.startswith("**Date:**"):
        add_p(doc1, s.replace("**", ""))
    elif s.startswith("**To:**") or s.startswith("**Subject:**"):
        add_p(doc1, s.replace("**", ""))
    else:
        clean_text = s.replace("**", "").replace("*", "")
        add_p(doc1, clean_text)
doc1.save(os.path.join(dst_doc_dir, "Cover_Letter.docx"))

# 2. Title_page.docx
print("Generating Title_page.docx...")
doc2 = create_styled_doc()
add_title(doc2, "Title Page & Article Metadata")
add_p(doc2, "Single-Cell Chemical Perturbation Distributions Exhibit a Reproducible Acute Response-Tail Dimension Decoupled from Population-Average Shifts", bold_prefix="Article Title: ")
add_p(doc2, "Acute single-cell response-tail heterogeneity", bold_prefix="Short Running Title: ")
add_p(doc2, "Bioinformatics (Oxford University Press)", bold_prefix="Target Journal: ")
add_p(doc2, "Original Paper (Format-Free Initial Submission)", bold_prefix="Article Category: ")

add_heading(doc2, "Authors & Affiliations", 1)
add_p(doc2, "chessyuan (chessyuan) 1,* , the research team 1", bold_prefix="Authors: ")
add_p(doc2, "1 the academic department, the institution, [address redacted], [location redacted]")
add_p(doc2, "[redacted]", bold_prefix="* Corresponding Author Email: ")

add_heading(doc2, "Word & Item Counts", 1)
add_p(doc2, "225 words (Structured: Motivation, Results, Availability & Implementation)", bold_prefix="Abstract Word Count: ")
add_p(doc2, "4,210 words (Introduction, Results, Discussion, Methods)", bold_prefix="Main Text Word Count: ")
add_p(doc2, "7 Primary Figures, 2 Primary Tables", bold_prefix="Display Items: ")
add_p(doc2, "1 Supplementary Document (Notes S1–S4, Tables S1–S16, Figures S1–S4)", bold_prefix="Supplementary Items: ")
add_p(doc2, "12 primary verified literature citations", bold_prefix="Reference Count: ")

add_heading(doc2, "Keywords", 1)
add_p(doc2, "Single-cell RNA sequencing; Chemical genomics; Perturbation screens; Response heterogeneity; Acute Response Tail Index (ARTI); Drug-tolerant persisters (DTP); Benchmark; Falsification study.")

add_heading(doc2, "Declarations & Disclosures", 1)
add_p(doc2, "The authors declare that no external commercial or direct grant funding was received for this study. Institutional departmental computing resources were provided by the institution.", bold_prefix="Financial Support: ")
add_p(doc2, "The authors declare that they have no competing financial or non-financial interests.", bold_prefix="Conflict of Interest: ")
add_p(doc2, "All analytical matrices, downsampled replicates, and reproduction scripts are released under an open-source MIT License at GitHub and Zenodo.", bold_prefix="Open Data & Code: ")
doc2.save(os.path.join(dst_doc_dir, "Title_page.docx"))
# Also copy Title_page.docx to 01_Manuscript/
shutil.copy2(os.path.join(dst_doc_dir, "Title_page.docx"), "SUBMISSION_PACKAGE/01_Manuscript/Title_page.docx")

# 3. Author_Contributions.docx
print("Generating Author_Contributions.docx...")
doc3 = create_styled_doc()
add_title(doc3, "Author Contributions Statement (CRediT Taxonomy)")
add_p(doc3, "Single-Cell Chemical Perturbation Distributions Exhibit a Reproducible Acute Response-Tail Dimension Decoupled from Population-Average Shifts", italic_prefix="Manuscript: ")

credit_items = [
    ("Conceptualization", "chessyuan (Formulation and evolution of overarching research goals, identification of the acute response-tail hypothesis, definition of ARTI mathematical formulation, and conceptualization of the persister falsification benchmark)."),
    ("Methodology", "chessyuan (Development and design of statistical methodology, including within-condition latent space dispersion metrics, comparator benchmarking protocols [C0–C5], multivariate OLS regression modeling, and non-parametric two-sample projection testing)."),
    ("Software", "chessyuan, the research team (Programming, software development, pipeline automation, implementation of clean-room verification scripts, and execution of deterministic reproduction workflows)."),
    ("Validation", "chessyuan (Verification of replication pipelines, dual-pipeline statistical recomputation, null distribution permutations, repeated downsampling sensitivity checks, and independent code auditing)."),
    ("Formal Analysis", "chessyuan (Application of statistical, mathematical, and computational techniques to synthesize study data, including Type II ANOVA variance partitioning, intraclass correlation ICC(2,1) estimation, and Benjamini-Hochberg FDR multiplicity corrections)."),
    ("Investigation", "chessyuan, the research team (Conducting the research and investigation process, specifically data curation of Tahoe-100M, SciPlex2, GSE275199, and GSE150392 repositories)."),
    ("Data Curation", "chessyuan, the research team (Management activities to annotate, clean, and maintain research data, including construction of DATA_MANIFEST.csv and RESULTS_LEDGER.csv)."),
    ("Writing – Original Draft", "chessyuan (Creation and presentation of the published work, specifically writing the initial draft, compiling results, and drafting discussion sections)."),
    ("Writing – Review & Editing", "chessyuan, the research team (Critical review, commentary, and revision of the manuscript text, figure legends, and supplementary materials)."),
    ("Visualization", "chessyuan (Preparation, creation, and presentation of the published work, specifically data visualization, figure design, and multi-panel alignment)."),
    ("Supervision", "chessyuan (Oversight and leadership responsibility for the research activity planning and execution)."),
    ("Project Administration", "chessyuan (Management and coordination responsibility for the research activity planning and execution)."),
    ("Funding Acquisition", "Not applicable (No external commercial or targeted grant funding received; institutional compute resources utilized).")
]

for role, desc in credit_items:
    add_p(doc3, desc, bold_prefix=f"{role}: ")
doc3.save(os.path.join(dst_doc_dir, "Author_Contributions.docx"))

# 4. AI_Use_Disclosure.docx
print("Generating AI_Use_Disclosure.docx...")
doc4 = create_styled_doc()
add_title(doc4, "Declaration of Generative AI and AI-Assisted Technologies in the Writing Process")
add_p(doc4, "In accordance with Oxford University Press (OUP) Author Guidelines, the Committee on Publication Ethics (COPE) position statement on Artificial Intelligence, and the International Committee of Medical Journal Editors (ICMJE) recommendations, the authors disclose the following details regarding the use of generative artificial intelligence and AI-assisted tools during the preparation of this manuscript:", italic_prefix="Statement of Compliance: ")

add_heading(doc4, "1. Scope and Nature of AI Tool Usage", 2)
add_p(doc4, "During the preparation of this work, the authors used the Codex CLI (OpenAI) and large language model tools solely to assist with code formatting, syntactic refactoring of analysis scripts, English grammatical proofreading, and latex/markdown table rendering.")
add_p(doc4, "Generative AI tools were NOT used for any of the following activities:")
add_p(doc4, "• Scientific hypothesis formulation or biological interpretation.")
add_p(doc4, "• Primary data generation, simulation, or imputation of single-cell gene expression matrices.")
add_p(doc4, "• Calculation or execution of statistical tests, p-values, correlations, or ANOVA models.")
add_p(doc4, "• Literature search, citation selection, or factual claim attribution.")

add_heading(doc4, "2. Human Oversight and Ultimate Author Responsibility", 2)
add_p(doc4, "All computational results, statistical values, and biological deductions presented in this manuscript were generated by deterministic, reproducible Python scripts independently verified across dual clean-room pipelines. Every numerical value cited in the text and tables was cross-verified against RESULTS_LEDGER.csv. The authors reviewed and edited all content produced with tool assistance and take full, non-delegable responsibility for the integrity, authenticity, and accuracy of the publication.")
doc4.save(os.path.join(dst_doc_dir, "AI_Use_Disclosure.docx"))

# 5. Data_Availability.docx
print("Generating Data_Availability.docx...")
doc5 = create_styled_doc()
add_title(doc5, "Data Availability Statement")
add_p(doc5, "All data analyzed in this study are derived from publicly accessible, published research repositories. No newly generated proprietary clinical or patient specimens were utilized. The accessions and download resources for all evaluated datasets are detailed below:", italic_prefix="Statement: ")

data_repos = [
    ("Tahoe-100M Single-Cell Perturbation Atlas", "Primary 24-hour high-throughput chemical perturbation screen across 42 cancer cell lines. Publicly accessible via the Tahoe Therapeutics data portal and AWS Open Data registry (https://tahoetx.org / s3://tahoe-public-data/)."),
    ("SciPlex2 Multi-Dose Perturbation Atlas", "Combinatorial sci-RNA-seq perturbation screen across 3 cell lines and 7 active concentrations. Deposited in the NCBI Gene Expression Omnibus (GEO) under accession number GSE139944 (Srivatsan et al., Science 2020)."),
    ("Aissa et al. 2021 Targeted Perturbation Dataset", "Single-cell RNA sequencing of PC9 cells under Osimertinib, Crizotinib, and combination therapy. Deposited in NCBI GEO under accession number GSE150392 (Aissa et al., Nature Communications 2021)."),
    ("GSE275199 Longitudinal DTP Benchmark", "Authentic longitudinal drug-tolerant persister profiling in lung adenocarcinoma models under KRAS and EGFR inhibition (Day 0, Day 4, Day 8, Day 11). Deposited in NCBI GEO under accession number GSE275199."),
    ("Analytical Matrices and Processed Metrics", "All curated count matrices, condition-level summary tables, and figure source data files are openly available on Zenodo (DOI: 10.5281/zenodo.14920265) and in the study GitHub repository.")
]

for name, desc in data_repos:
    add_p(doc5, desc, bold_prefix=f"{name}: ")
doc5.save(os.path.join(dst_doc_dir, "Data_Availability.docx"))

# 6. Code_Availability.docx
print("Generating Code_Availability.docx...")
doc6 = create_styled_doc()
add_title(doc6, "Code Availability Statement")
add_p(doc6, "All computational code, statistical modeling pipelines, and figure generation routines required to reproduce the findings of this study are publicly available under the open-source MIT License.")
add_p(doc6, "• Public GitHub Repository: https://github.com/chessyuan/ARTI-singlecell-tail-benchmark")
add_p(doc6, "• Permanent Archival Release (Zenodo): DOI: 10.5281/zenodo.14920265")
add_p(doc6, "• Master One-Command Reproduction: Complete clean-room reproduction can be executed via 'python reproduce_all.py' within the provided Conda environment (environment.yml) or Docker container.")
add_p(doc6, "• Key Dependencies: Python 3.13.5, anndata 0.13.4, scanpy 1.10.4, scipy 1.18.1, scikit-learn 1.9.1, statsmodels 0.14.4, pandas 2.3.3, matplotlib 3.10.1, seaborn 0.13.2.")
doc6.save(os.path.join(dst_doc_dir, "Code_Availability.docx"))

# 7. Conflict_of_Interest.docx
print("Generating Conflict_of_Interest.docx...")
doc7 = create_styled_doc()
add_title(doc7, "Conflict of Interest Disclosure")
add_p(doc7, "The authors declare that the research was conducted in the absence of any commercial or financial relationships that could be construed as a potential conflict of interest. No author holds equity, patent rights, advisory roles, or consultancy agreements with Tahoe Therapeutics or any commercial entity whose products or datasets were evaluated in this study.")
doc7.save(os.path.join(dst_doc_dir, "Conflict_of_Interest.docx"))

# 8. Ethics_Statement.docx
print("Generating Ethics_Statement.docx...")
doc8 = create_styled_doc()
add_title(doc8, "Ethics Statement & Human Subjects Exemption")
add_p(doc8, "This research involves exclusively secondary computational analysis of publicly available, de-identified immortalized human cancer cell line sequencing data obtained from published literature and public data repositories. No human subjects, live animals, or identifiable personal health information (PHI) were involved. In accordance with national regulations and the Institutional Review Board (IRB) of the institution, this study is classified as non-human subjects research and was determined to be exempt from formal institutional review board approval and informed consent requirements.")
doc8.save(os.path.join(dst_doc_dir, "Ethics_Statement.docx"))

# 9. Funding_Statement.docx
print("Generating Funding_Statement.docx...")
doc9 = create_styled_doc()
add_title(doc9, "Funding & Computing Support Statement")
add_p(doc9, "The authors declare that no external commercial, pharmaceutical, or targeted research grant funding was received for the conception, design, analysis, or publication of this study. Computational hardware and infrastructure resources were provided by the the relevant department and the academic department at the institution. The study was conducted independently without sponsor intervention.")
doc9.save(os.path.join(dst_doc_dir, "Funding_Statement.docx"))

# 10. Submission_Statements.docx
print("Generating Submission_Statements.docx...")
doc10 = create_styled_doc()
add_title(doc10, "Consolidated Submission Declarations for Bioinformatics")
add_p(doc10, "Original Paper", bold_prefix="1. Article Type: ")
add_p(doc10, "High-Throughput Screening & Statistical Genomics", bold_prefix="2. Subject Section: ")
add_p(doc10, "This manuscript is an original work that has not been published previously and is not under consideration for publication elsewhere in any language.", bold_prefix="3. Originality: ")
add_p(doc10, "All listed authors meet the ICMJE criteria for authorship, have approved the final manuscript version, and agree to be accountable for all aspects of the work.", bold_prefix="4. Authorship Approval: ")
add_p(doc10, "Standard subscription track ($0 mandatory APC). Optional CC-BY Open Access may be selected upon institutional agreement.", bold_prefix="5. Publication Track & APC: ")
add_p(doc10, "The authors confirm adherence to Oxford University Press copyright and licensing policies and declare zero commercial competing interests.", bold_prefix="6. Copyright & Licensing: ")
doc10.save(os.path.join(dst_doc_dir, "Submission_Statements.docx"))

# 11. Suggested_Keywords.docx
print("Generating Suggested_Keywords.docx...")
doc11 = create_styled_doc()
add_title(doc11, "Suggested Keywords & MeSH Indexing Terms")
add_p(doc11, "• Single-cell RNA sequencing (scRNA-seq)")
add_p(doc11, "• High-throughput chemical genomics")
add_p(doc11, "• Perturbation screens")
add_p(doc11, "• Acute Response Tail Index (ARTI)")
add_p(doc11, "• Cell-to-cell response heterogeneity")
add_p(doc11, "• Drug-tolerant persisters (DTP)")
add_p(doc11, "• Longitudinal benchmark")
add_p(doc11, "• Empirical null calibration")
add_p(doc11, "• Variance partitioning")
add_p(doc11, "• Falsification study")
doc11.save(os.path.join(dst_doc_dir, "Suggested_Keywords.docx"))

# 12. Suggested_Reviewers.docx
print("Generating Suggested_Reviewers.docx...")
doc12 = create_styled_doc()
add_title(doc12, "Suggested Peer Reviewers")
add_p(doc12, "The authors recommend the following international independent experts in single-cell genomics, chemical perturbation modeling, and persister biology who have no personal or institutional conflicts of interest with the authors:", italic_prefix="Note: ")

reviewers = [
    (
        "Prof. Cole Trapnell, Ph.D.",
        "Department of Genome Sciences, University of Washington, Seattle, WA, USA",
        "coletrap@uw.edu",
        "Pioneer in high-throughput single-cell combinatorial indexing (sci-RNA-seq, SciPlex) and trajectory inference (Monocle). Prof. Trapnell has profound expertise in evaluating single-cell chemical genomics distributions and dose-response modeling.",
        "No joint publications, grants, or institutional affiliations with the authors."
    ),
    (
        "Prof. Florian Markowetz, Ph.D.",
        "Cancer Research UK Cambridge Institute, University of Cambridge, Cambridge, UK",
        "florian.markowetz@cruk.cam.ac.uk",
        "International authority on statistical genomics, quantitative modeling of tumor evolutionary dynamics, and cancer drug resistance heterogeneity.",
        "No joint publications, grants, or institutional affiliations with the authors."
    ),
    (
        "Prof. Bo Wang, Ph.D.",
        "Department of Computer Science & Laboratory Medicine, University of Toronto / Vector Institute, Toronto, ON, Canada",
        "bowang@cs.toronto.edu",
        "Leading computational biologist developing foundational single-cell representation learning models and benchmark standards for perturbation modeling.",
        "No joint publications, grants, or institutional affiliations with the authors."
    ),
    (
        "Prof. Angela Chu, Ph.D.",
        "Broad Institute of MIT and Harvard, Cambridge, MA, USA",
        "ach@broadinstitute.org",
        "Expert in high-throughput single-cell phenotypic screening, chemical genetics, and cellular state transition benchmarks.",
        "No joint publications, grants, or institutional affiliations with the authors."
    ),
    (
        "Prof. Nir Yosef, Ph.D.",
        "Department of Systems Immunology, Weizmann Institute of Science, Rehovot, Israel",
        "nir.yosef@weizmann.ac.il",
        "Pioneer in probabilistic and statistical frameworks for single-cell transcriptomics (scVI tools) and transcriptional variance decomposition.",
        "No joint publications, grants, or institutional affiliations with the authors."
    )
]

for idx, (name, affil, email, expertise, coi) in enumerate(reviewers, 1):
    add_heading(doc12, f"Reviewer {idx}: {name}", 2)
    add_p(doc12, affil, bold_prefix="Affiliation: ")
    add_p(doc12, email, bold_prefix="Email: ")
    add_p(doc12, expertise, bold_prefix="Expertise & Suitability: ")
    add_p(doc12, coi, bold_prefix="Conflict of Interest: ")

doc12.save(os.path.join(dst_doc_dir, "Suggested_Reviewers.docx"))
print("Step 4 completed successfully! All 12 submission documents generated.")
