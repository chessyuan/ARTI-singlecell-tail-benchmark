import urllib.request
import os

base_url = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12063144/bin/'
files = ['mmc1.pdf', 'mmc2.xlsx', 'mmc3.xlsx', 'mmc4.xlsx']

for f in files:
    out_path = os.path.join('data_raw/external_dtp/gse275199', f)
    if not os.path.exists(out_path):
        print(f'Downloading {f}...')
        req = urllib.request.Request(base_url + f, headers={'User-Agent': 'Mozilla/5.0'})
        try:
            with urllib.request.urlopen(req) as resp, open(out_path, 'wb') as out_f:
                out_f.write(resp.read())
            print(f'Done {f}, size = {os.path.getsize(out_path)} bytes')
        except Exception as e:
            print(f'Failed {f}: {e}')
    else:
        print(f'Already exists: {f}')
