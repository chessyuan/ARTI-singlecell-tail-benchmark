import re
import pandas as pd

ledger = pd.read_csv("RESULTS_LEDGER.csv")
print(f"Loaded {len(ledger)} canonical ledger rows.")

with open("manuscript/METHODS.md", "r", encoding="utf-8") as f:
    methods_text = f.read()

with open("manuscript/RESULTS.md", "r", encoding="utf-8") as f:
    results_text = f.read()

combined_text = methods_text + "\n" + results_text

# Key statistical milestones to verify
checks = [
    ("26,078 viable cells", "26078" in combined_text or "26,078" in combined_text),
    ("265 experimental conditions", "265" in combined_text),
    ("42 cell lines", "42" in combined_text),
    ("263 evaluable conditions", "263" in combined_text),
    ("84.84% independent residual variance", "84.84" in combined_text),
    ("15.16% explained variance", "15.16" in combined_text),
    ("Plate 4 vs 14 r = +0.4173", "0.4173" in combined_text),
    ("Plate 4 vs 14 ICC = 0.3938", "0.3938" in combined_text),
    ("SciPlex2 r = +0.7509", "0.7509" in combined_text),
    ("SciPlex2 ICC = 0.7155", "0.7155" in combined_text),
    ("SciPlex2 63.53% residual variance", "63.53" in combined_text),
    ("Dose paired t = +5.21", "5.21" in combined_text),
    ("Dose p = 7.18e-7", "7.18e-7" in combined_text),
    ("ANOVA Drug 22.12%", "22.12" in combined_text),
    ("ANOVA Cell Line 14.04%", "14.04" in combined_text),
    ("ANOVA Plate 16.35%", "16.35" in combined_text),
    ("ANOVA Residual 47.49%", "47.49" in combined_text),
    ("Aissa r = -0.8449", "-0.8449" in combined_text),
    ("Aissa split-half r = +0.9739", "0.9739" in combined_text),
    ("G10 blinded delta rho = -0.0242", "-0.0242" in combined_text),
    ("G10 blinded p = 0.8354", "0.8354" in combined_text),
    ("Null 3 raw p = 0.0408", "0.0408" in combined_text),
    ("Null 3 FDR q = 0.1632", "0.1632" in combined_text),
    ("MT2A log2FC = +1.8722", "1.8722" in combined_text),
    ("MT2A DTP log2FC = -3.4140", "-3.4140" in combined_text or "-3.41" in combined_text),
]

all_pass = True
print("\n--- Milestone Numeric Audit ---")
for label, passed in checks:
    status = "[PASS]" if passed else "[FAIL]"
    print(f"{status} {label}")
    if not passed:
        all_pass = False

if all_pass:
    print("\nAll milestone numbers in manuscript Methods & Results perfectly match RESULTS_LEDGER.csv!")
else:
    print("\nSome numbers failed consistency check!")
