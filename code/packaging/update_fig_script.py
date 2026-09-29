import os
import shutil
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pypdf import PdfWriter, PdfReader

os.makedirs("SUBMISSION_PACKAGE/02_Figures", exist_ok=True)
os.makedirs("figures/main", exist_ok=True)

# 1. Update code/figures/generate_all_main_figures.py to save both png and pdf
with open("code/figures/generate_all_main_figures.py", "r", encoding="utf-8") as f:
    code = f.read()

# Replace png-only saves with both png and pdf
for i in range(1, 8):
    old_save = f"fig{i}_"
    # Find matching savefig
    # We will do replacement for each figure save
code_mod = code
fig_names = [
    (1, "fig1_art_framework_qc_null"),
    (2, "fig2_arti_vs_average_response"),
    (3, "fig3_replicate_dose_reproducibility"),
    (4, "fig4_pan_tahoe_landscape_variance"),
    (5, "fig5_external_acute_validation"),
    (6, "fig6_metric_benchmark_competition"),
    (7, "fig7_dtp_negative_boundary")
]

for num, name in fig_names:
    old_str = f"plt.savefig('figures/main/{name}.png', dpi=300)"
    new_str = f"plt.savefig('figures/main/{name}.png', dpi=300)\nplt.savefig('figures/main/{name}.pdf')\nplt.savefig('figures/main/{name}.tif', dpi=300)"
    if old_str in code_mod:
        code_mod = code_mod.replace(old_str, new_str)
    else:
        print(f"Warning: {old_str} not found in code")

with open("code/figures/generate_all_main_figures.py", "w", encoding="utf-8") as f:
    f.write(code_mod)

print("Updated code/figures/generate_all_main_figures.py with PDF and TIF exports.")
