import re

with open('data_raw/external_dtp/gse275199/PMC12063144.xml', 'r', encoding='utf-8') as f:
    text = f.read()

matches = re.findall(r'(https?://[^\s"\'<>]+)', text)
repos = [m for m in matches if any(x in m.lower() for x in ['github', 'zenodo', 'figshare', 'dryad', 'osf', 'synapse', 'geo', 'ncbi'])]
print('Relevant links found:')
for r in set(repos):
    print(' ', r)

avail = re.findall(r'<sec[^>]*id="sec[0-9.]*"[^>]*>(.*?)</sec>', text, re.DOTALL)
for a in avail:
    if 'availability' in a.lower():
        print('=== Section ===')
        print(re.sub(r'<[^>]+>', ' ', a).strip()[:500])
