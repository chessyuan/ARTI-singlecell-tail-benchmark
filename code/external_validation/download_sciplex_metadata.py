import os
import urllib.request

out_dir = 'data_raw/external_acute/sciplex2'
os.makedirs(out_dir, exist_ok=True)

base_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139944/suppl/'
files_to_download = [
    'GSM4150377_sciPlex2_hashSampleSheet.txt.gz',
    'GSM4150377_sciPlex2_pData.txt.gz',
    'GSM4150377_sciPlex2_A549_Transcription_Modulators_cell.annotations.txt.gz',
    'GSM4150377_sciPlex2_A549_Transcription_Modulators_gene.annotations.txt.gz'
]

for fname in files_to_download:
    dest = os.path.join(out_dir, fname)
    if not os.path.exists(dest):
        url = base_url + fname
        print(f"Downloading {fname}...")
        urllib.request.urlretrieve(url, dest)
        print(f"Downloaded {fname}: {os.path.getsize(dest)} bytes.")
    else:
        print(f"Already exists: {fname}")
