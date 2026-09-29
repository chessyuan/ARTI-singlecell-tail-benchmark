import re

with open('data_raw/external_dtp/gse275199/PMC12063144.xml', 'r', encoding='utf-8') as f:
    text = f.read()

supps = re.findall(r'<supplementary-material[^>]*>(.*?)</supplementary-material>', text, re.DOTALL)
print('Number of supplementary materials:', len(supps))
for s in supps:
    clean = re.sub(r'<[^>]+>', ' ', s)
    hrefs = re.findall(r'xlink:href=[\'"]([^\'"]+)[\'"]', s)
    print('  supp:', hrefs, clean[:200].strip())
