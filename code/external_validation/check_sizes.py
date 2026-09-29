import urllib.request

base_url = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4150nnn/GSM4150377/suppl/'
files = [
    'GSM4150377_sciPlex2_A549_Transcription_Modulators_cell.annotations.txt.gz',
    'GSM4150377_sciPlex2_A549_Transcription_Modulators_gene.annotations.txt.gz',
    'GSM4150377_sciPlex2_hashSampleSheet.txt.gz',
    'GSM4150377_sciPlex2_pData.txt.gz',
    'GSM4150377_sciPlex2_A549_Transcription_Modulators_UMI.count.matrix.gz'
]

for f in files:
    url = base_url + f
    req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            cl = resp.headers.get('Content-Length')
            mb = int(cl)/(1024*1024) if cl else 0
            print(f"{f}: {cl} bytes ({mb:.2f} MB)")
    except Exception as e:
        print(f"{f}: error {e}")
