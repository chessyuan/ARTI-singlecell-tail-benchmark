import os
import shutil
import docx
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
import win32com.client

dst_fig_dir = os.path.abspath("SUBMISSION_PACKAGE/02_Figures")
os.makedirs(dst_fig_dir, exist_ok=True)

fig_mapping = [
    (1, "fig1_art_framework_qc_null", "Figure_1"),
    (2, "fig2_arti_vs_average_response", "Figure_2"),
    (3, "fig3_replicate_dose_reproducibility", "Figure_3"),
    (4, "fig4_pan_tahoe_landscape_variance", "Figure_4"),
    (5, "fig5_external_acute_validation", "Figure_5"),
    (6, "fig6_metric_benchmark_competition", "Figure_6"),
    (7, "fig7_dtp_negative_boundary", "Figure_7")
]

# Copy individual files
for num, src_name, dst_name in fig_mapping:
    for ext in [".png", ".pdf", ".tif"]:
        src_path = os.path.join("figures/main", src_name + ext)
        dst_path = os.path.join(dst_fig_dir, dst_name + ext)
        shutil.copy2(src_path, dst_path)
        print(f"Copied {src_path} -> {dst_path} ({os.path.getsize(dst_path):,} bytes)")

# Parse figure legends from MAIN_MANUSCRIPT.md
with open("submission/final/MAIN_MANUSCRIPT.md", "r", encoding="utf-8") as f:
    ms_text = f.read()

pos_fig = ms_text.find("## Figure Legends")
pos_supp = ms_text.find("## Supplementary Information")
legends_block = ms_text[pos_fig:pos_supp].strip()

# Split legends by Figure X
import re
fig_chunks = re.split(r'\n(?=\*\*Figure \d+\.)', legends_block)
legends = {}
for chunk in fig_chunks:
    m = re.match(r'\*\*Figure (\d+)\.\s*(.*?)\*\*(.*)', chunk, re.DOTALL)
    if m:
        fig_idx = int(m.group(1))
        title = m.group(2).strip()
        body = m.group(3).strip()
        legends[fig_idx] = (title, body)

print(f"Parsed legends for figures: {list(legends.keys())}")

# Create Word document for Figures_combined_reference
doc = docx.Document()
for s in doc.sections:
    s.top_margin = Inches(0.8)
    s.bottom_margin = Inches(0.8)
    s.left_margin = Inches(0.8)
    s.right_margin = Inches(0.8)

title_p = doc.add_paragraph()
title_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r_t = title_p.add_run("Combined Primary Figures and Legends")
r_t.font.name = "Times New Roman"
r_t.font.size = Pt(16)
r_t.font.bold = True

sub_p = doc.add_paragraph()
sub_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
r_s = sub_p.add_run("Bioinformatics — Original Paper Submission\nSingle-Cell Chemical Perturbation Distributions Exhibit a Reproducible Acute Response-Tail Dimension Decoupled from Population-Average Shifts\nchessyuan et al., the institution\n")
r_s.font.name = "Times New Roman"
r_s.font.size = Pt(11)
r_s.font.italic = True

for num, src_name, dst_name in fig_mapping:
    doc.add_page_break()
    png_path = os.path.abspath(os.path.join(dst_fig_dir, f"{dst_name}.png"))
    
    # Add figure heading
    title, body = legends.get(num, (f"Figure {num}", ""))
    hp = doc.add_paragraph()
    hp.paragraph_format.space_before = Pt(0)
    hp.paragraph_format.space_after = Pt(6)
    r_h = hp.add_run(f"Figure {num}. {title}")
    r_h.font.name = "Times New Roman"
    r_h.font.size = Pt(12)
    r_h.font.bold = True
    
    # Add image
    img_p = doc.add_paragraph()
    img_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    img_p.paragraph_format.space_before = Pt(4)
    img_p.paragraph_format.space_after = Pt(8)
    img_run = img_p.add_run()
    # Choose size to fit page cleanly
    img_run.add_picture(png_path, width=Inches(6.8))
    
    # Add legend body
    leg_p = doc.add_paragraph()
    leg_p.paragraph_format.space_before = Pt(4)
    leg_p.paragraph_format.space_after = Pt(12)
    leg_p.paragraph_format.line_spacing = 1.15
    r_body = leg_p.add_run(body)
    r_body.font.name = "Times New Roman"
    r_body.font.size = Pt(9.5)

combo_docx = os.path.abspath("SUBMISSION_PACKAGE/02_Figures/Figures_combined_reference.docx")
combo_pdf = os.path.abspath("SUBMISSION_PACKAGE/02_Figures/Figures_combined_reference.pdf")
doc.save(combo_docx)
print(f"Saved {combo_docx}")

# Convert to PDF via Word COM
word = win32com.client.DispatchEx("Word.Application")
word.Visible = False
try:
    doc_com = word.Documents.Open(combo_docx)
    doc_com.SaveAs(combo_pdf, FileFormat=17) # 17 = wdFormatPDF
    doc_com.Close()
    print(f"Exported combined PDF via Word COM: {combo_pdf} ({os.path.getsize(combo_pdf):,} bytes)")
finally:
    word.Quit()
