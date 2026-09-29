import os
import sys
import glob
import re

# Core forbidden patterns (unconditional)
FORBIDDEN_PATTERNS = [
    r"first\s+discovery\s+of\s+heterogeneity",
    r"first\s+demonstration\s+that\s+single-cell\s+drug\s+responses\s+are\s+heterogeneous",
    r"mathematically\s+orthogonal",
    r"independent\s+biological\s+information",
    r"biological\s+residual\s+variance",
    r"84\.84%\s+biological\s+variance",
    r"84\.84%\s+independent\s+information",
    r"proves?\s+(the\s+residual\s+is\s+)?biological",
    r"proves?\s+biological\s+signal",
    r"predicts?\s+persisters?",
    r"predicts?\s+clinical\s+resistance",
    r"predicts?\s+DTP",
    r"predicts?\s+future\s+survival",
    r"adaptive\s+survival\s+tail",
    r"adaptive\s+survival\s+state",
    r"represents?\s+survival",
    r"cytotoxicity",
    r"cell\s+killing",
    r"Adaptive\s+Heterogeneity\s+Window",
    r"DTP\s+precursor",
    r"synergy\b",
    r"cures?\s+DTP",
    r"combination\s+cures",
    r"persisters?\s+discovered\s+in\s+Tahoe",
    r"genetic\s+background",
    r"drug\s+pharmacology",
]

# Context-sensitive patterns
# 'viable cells' or 'surviving cells' are only allowed if explicitly discussing the lack/absence of viability assays
CONTEXT_SENSITIVE_PATTERNS = [
    (r"\bviable\s+cells?\b", [r"lack\s+of\s+functional\s+viability", r"without\s+functional\s+viability", r"no\s+direct\s+viability\s+assay", r"viability\s+assays?"]),
    (r"\bsurviving\s+cells?\b", [r"lack\s+of\s+functional\s+viability"]),
    (r"\bliving\s+cells?\b", [r"lack\s+of\s+functional\s+viability"]),
]

def scan_markdown_files(directory="manuscript"):
    md_files = glob.glob(os.path.join(directory, "*.md"))
    if not md_files:
        print(f"[*] No markdown files found in {directory}/ to scan yet.")
        return True
        
    violations = []
    
    for filepath in sorted(md_files):
        with open(filepath, "r", encoding="utf-8") as f:
            lines = f.readlines()
            
        in_code_block = False
        for line_num, line in enumerate(lines, 1):
            if line.strip().startswith("```"):
                in_code_block = not in_code_block
                continue
            if in_code_block:
                continue
                
            # Skip lines in markdown that explicitly discuss forbidden rules, negative boundaries, or audit notes
            if any(marker in line for marker in ["FORBIDDEN", "Forbidden Claims", "NEVER termed", "Negative Boundary", "cautioning against interpreting", "prohibited", "prohibition", "Constraints"]) and not ("Result" in line):
                continue
                
            for pattern in FORBIDDEN_PATTERNS:
                if re.search(pattern, line, re.IGNORECASE):
                    violations.append((filepath, line_num, line.strip(), pattern))
                    
            for pattern, allowed_exceptions in CONTEXT_SENSITIVE_PATTERNS:
                if re.search(pattern, line, re.IGNORECASE):
                    if not any(re.search(exc, line, re.IGNORECASE) for exc in allowed_exceptions):
                        violations.append((filepath, line_num, line.strip(), f"Context-sensitive: {pattern}"))
                    
    if violations:
        print(f"\n[FAIL] Found {len(violations)} forbidden claim violations in {directory}/:")
        for fp, ln, text, pat in violations:
            print(f"  {fp}:{ln} -> Matches pattern: '{pat}'")
            print(f"    Line text: \"{text[:120]}\"")
        return False
    else:
        print(f"\n[PASS] Scanned {len(md_files)} markdown files in {directory}/. Zero forbidden claim violations detected.")
        return True

if __name__ == "__main__":
    success = scan_markdown_files()
    sys.exit(0 if success else 1)
