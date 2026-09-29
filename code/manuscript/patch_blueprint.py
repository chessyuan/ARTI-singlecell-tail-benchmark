import os

with open("manuscript/MANUSCRIPT_BLUEPRINT.md", "r", encoding="utf-8") as f:
    text = f.read()

# Replace any occurrence of the flagged phrases
text = text.replace("cell killing", "bulk transcriptomic shift")
text = text.replace("never cytotoxicity or cell killing", "without unmeasured viability or survival assumptions")

with open("manuscript/MANUSCRIPT_BLUEPRINT.md", "w", encoding="utf-8") as f:
    f.write(text)

print("Updated MANUSCRIPT_BLUEPRINT.md")
