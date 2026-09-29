import urllib.request
import re

for gsm in ['GSM8472663', 'GSM8472664']:
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req) as resp:
            html = resp.read().decode('utf-8')
            print(f'=== {gsm} ===')
            # Extract table rows
            rows = re.findall(r'<tr valign="top"><td[^>]*>(.*?)</td><td[^>]*>(.*?)</td></tr>', html, re.DOTALL)
            for k, v in rows:
                k_clean = re.sub(r'<[^>]+>', '', k).strip()
                v_clean = re.sub(r'<[^>]+>', '', v).strip()
                if any(x in k_clean.lower() for x in ['title', 'source', 'organism', 'characteristics', 'growth', 'molecule', 'extract', 'description', 'supplementary']):
                    print(f'  {k_clean}: {v_clean[:120]}')
    except Exception as e:
        print(f'{gsm} error:', e)
