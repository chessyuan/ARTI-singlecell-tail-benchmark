import urllib.request
import re

url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE275nnn/GSE275199/'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as resp:
        html = resp.read().decode('utf-8')
        links = re.findall(r'href=[\'\"]([^\'\"]+)[\'\"]', html)
        print('Files in GSE275199 root:', links)
except Exception as e:
    print('Error:', e)
