import os
import shutil
import subprocess
import docx
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import parse_xml, OxmlElement
from docx.oxml.ns import nsdecls, qn
import win32com.client

dst_ms_dir = os.path.abspath("SUBMISSION_PACKAGE/01_Manuscript")
os.makedirs(dst_ms_dir, exist_ok=True)

src_md = os.path.abspath("submission/final/MAIN_MANUSCRIPT.md")
temp_docx = os.path.abspath("temp_pandoc_ms.docx")

# Step 1: Run Pandoc to convert Markdown to DOCX with native typography & math
cmd = ["pandoc", src_md, "-o", temp_docx]
print("Running pandoc...", cmd)
subprocess.run(cmd, check=True)

# Step 2: Post-process with python-docx
print("Post-processing Manuscript_full.docx...")
doc = docx.Document(temp_docx)

# Configure margins, headers, footers, line numbers for all sections
for section in doc.sections:
    section.top_margin = Inches(1.0)
    section.bottom_margin = Inches(1.0)
    section.left_margin = Inches(1.0)
    section.right_margin = Inches(1.0)
    
    # Add continuous line numbering
    sectPr = section._sectPr
    # Remove any existing lnNumType if present
    for existing in sectPr.findall(qn('w:lnNumType')):
        sectPr.remove(existing)
    xml_str = f'<w:lnNumType {nsdecls("w")} w:countBy="1" w:restart="continuous"/>'
    sectPr.append(parse_xml(xml_str))
    
    # Header: Running title
    header = section.header
    hp = header.paragraphs[0]
    hp.text = "Bioinformatics | Acute single-cell response-tail heterogeneity"
    hp.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    if hp.runs:
        hp.runs[0].font.name = "Times New Roman"
        hp.runs[0].font.size = Pt(8.5)
        hp.runs[0].font.italic = True
        hp.runs[0].font.color.rgb = RGBColor(128, 128, 128)
        
    # Footer: Page number
    footer = section.footer
    fp = footer.paragraphs[0]
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    frun = fp.add_run("Page ")
    frun.font.name = "Times New Roman"
    frun.font.size = Pt(9)
    # Add page field
    fldSimple = parse_xml(r'<w:fldSimple %s w:instr="PAGE"/>' % nsdecls('w'))
    fp._p.append(fldSimple)

# Style paragraphs
for p in doc.paragraphs:
    # Set line spacing and font
    if p.style.name.startswith("Heading 1"):
        p.paragraph_format.space_before = Pt(14)
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.keep_with_next = True
        for r in p.runs:
            r.font.name = "Times New Roman"
            r.font.size = Pt(14)
            r.font.bold = True
            r.font.color.rgb = RGBColor(0, 0, 0)
    elif p.style.name.startswith("Heading 2"):
        p.paragraph_format.space_before = Pt(12)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.keep_with_next = True
        for r in p.runs:
            r.font.name = "Times New Roman"
            r.font.size = Pt(12)
            r.font.bold = True
            r.font.color.rgb = RGBColor(0, 0, 0)
    elif p.style.name.startswith("Heading 3"):
        p.paragraph_format.space_before = Pt(10)
        p.paragraph_format.space_after = Pt(3)
        p.paragraph_format.keep_with_next = True
        for r in p.runs:
            r.font.name = "Times New Roman"
            r.font.size = Pt(11)
            r.font.bold = True
            r.font.color.rgb = RGBColor(0, 0, 0)
    else:
        # Body text
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.line_spacing = 1.35
        for r in p.runs:
            r.font.name = "Times New Roman"
            if not r.font.size:
                r.font.size = Pt(11)

full_docx = os.path.join(dst_ms_dir, "Manuscript_full.docx")
doc.save(full_docx)
print(f"Saved {full_docx} ({os.path.getsize(full_docx):,} bytes)")

# Step 3: Create Blinded Manuscript
print("Creating Manuscript_blinded.docx...")
doc_blind = docx.Document(full_docx)
# Anonymize author info in first 10 paragraphs
for p in doc_blind.paragraphs[:15]:
    t = p.text
    if "chessyuan" in t or "the institution" in t or "[redacted]" in t or "Working Group" in t:
        p.text = "[Author details, institutional affiliations, and contact information removed for double-blind peer review]"
        p.runs[0].font.name = "Times New Roman"
        p.runs[0].font.size = Pt(11)
        p.runs[0].font.italic = True

blind_docx = os.path.join(dst_ms_dir, "Manuscript_blinded.docx")
doc_blind.save(blind_docx)
print(f"Saved {blind_docx} ({os.path.getsize(blind_docx):,} bytes)")

# Step 4: Export Reference PDF via Word COM
print("Exporting Manuscript_reference.pdf via Word COM...")
ref_pdf = os.path.join(dst_ms_dir, "Manuscript_reference.pdf")
word = win32com.client.DispatchEx("Word.Application")
word.Visible = False
try:
    d = word.Documents.Open(full_docx)
    d.SaveAs(ref_pdf, FileFormat=17) # 17 = wdFormatPDF
    d.Close()
    print(f"Exported Manuscript_reference.pdf ({os.path.getsize(ref_pdf):,} bytes)")
finally:
    word.Quit()

if os.path.exists(temp_docx):
    os.remove(temp_docx)

print("Step 5 completed successfully!")
