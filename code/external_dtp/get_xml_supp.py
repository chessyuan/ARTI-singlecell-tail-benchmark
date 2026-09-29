import re

with open('data_raw/external_dtp/gse275199/PMC12063144.xml', 'r', encoding='utf-8') as f:
    text = f.read()

# Find all occurrences of mmc
matches = re.findall(r'.{0,50}mmc[0-9].{0,100}', text)
for m in set(matches):
    print(m)
