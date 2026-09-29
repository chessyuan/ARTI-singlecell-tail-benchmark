with open("STATUS.md", "r", encoding="utf-8") as f:
    status_content = f.read()

status_update = """
CURRENT_PHASE:
PHASE M (MANUSCRIPT ENGINEERING) - M1-M4 COMPLETE

ANALYSIS_FREEZE_STATUS:
LOCKED (Git commit e226c1e998ca079244bf9dfbd3c1f0535d6aa6c7, Tag v1.0-analysis-frozen; audit/ANALYSIS_FREEZE.md)

MANUSCRIPT_BLUEPRINT:
ACTIVE (manuscript/MANUSCRIPT_BLUEPRINT.md, FIGURE_RESULT_MAP.csv, TABLE_RESULT_MAP.csv, CLAIM_SENTENCE_MAP.csv)

METHODS_STATUS:
COMPLETE (manuscript/METHODS.md, 2,543 words, 12 sections, 100% verified against RESULTS_LEDGER.csv)

RESULTS_STATUS:
COMPLETE (manuscript/RESULTS.md, 2,343 words, 7 modules matching Figures 1-7, 100% verified)

FORBIDDEN_CLAIMS_SCAN:
PASS (0 violations detected by code/manuscript/check_forbidden_claims.py)

NUMERIC_AUDIT_STATUS:
PASS (100% concordance verified by code/manuscript/verify_m1_m4_consistency.py)

POST_FREEZE_AMENDMENT_REQUIRED:
NO (Zero bugs, zero source data inconsistencies, zero statistical unit errors)

NEXT_ACTION:
Phase M5 (Introduction, Discussion, Abstract, Title candidate formulation)

LAST_UPDATED:
2026-09-29
"""

status_content = status_content.replace("GO (MANUSCRIPT_READY: All 28 P13 sub-sections audited; 24 questions resolved; clean-room reproduction verified)", "GO (PHASE_M1_M4_COMPLETE: Analysis frozen, Blueprint locked, Methods and Results drafted)")
status_content = status_content.replace("Draft manuscript following strict Allowed Claims policy of audit/FINAL_CLAIM_LOCK.md.", "Execute Phase M5 (Introduction, Discussion, Abstract, Titles) following manuscript blueprint.")

with open("STATUS.md", "w", encoding="utf-8") as f:
    f.write(status_content.strip() + "\n" + status_update.strip() + "\n")

with open("DECISIONS.md", "r", encoding="utf-8") as f:
    dec_content = f.read()

dec_update = """
- **DEC-019 (2026-09-29)**: Phase M1–M4 Manuscript Engineering Milestone Complete:
  - Formally locked Analysis Freeze in `audit/ANALYSIS_FREEZE.md` (Git commit `e226c1e`, Git tag `v1.0-analysis-frozen`).
  - Established `manuscript/MANUSCRIPT_BLUEPRINT.md` defining central thesis, story architecture, and 7 result modules.
  - Generated complete mapping matrices: `manuscript/FIGURE_RESULT_MAP.csv` (7 figures, 25 panels mapped to result IDs), `manuscript/TABLE_RESULT_MAP.csv` (16 supplementary tables mapped to results), and `manuscript/CLAIM_SENTENCE_MAP.csv` (17 core claims, all ALLOWED or QUALIFIED).
  - Deployed automated forbidden claim scanner (`code/manuscript/check_forbidden_claims.py`); zero violations detected across all manuscript drafts.
  - Drafted comprehensive `manuscript/METHODS.md` (2,543 words across 12 mathematical and computational sections).
  - Drafted complete `manuscript/RESULTS.md` (2,343 words across 7 result modules matching Figures 1–7).
  - Verified 100% numerical consistency between manuscript text and `RESULTS_LEDGER.csv` via automated audit script (`code/manuscript/verify_m1_m4_consistency.py`).
  - Confirmed 0 numeric inconsistencies, 0 forbidden claims, and no requirement for `POST_FREEZE_AMENDMENT.md`. Manuscript readiness stands at Phase M5.
"""

with open("DECISIONS.md", "w", encoding="utf-8") as f:
    f.write(dec_content.strip() + "\n" + dec_update.strip() + "\n")

print("Updated STATUS.md and DECISIONS.md")
