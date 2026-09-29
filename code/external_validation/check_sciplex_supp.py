import urllib.request
import re

ftp_url = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4150nnn/GSM4150377/suppl/'
req = urllib.request.Request(ftp_url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req, timeout=15) as resp:
        html = resp.read().decode('utf-8')
        links = re.findall(r'href="([^"]+)"', html)
        print("GSM4150377 Supplementary files:")
        for l in links:
            if not l.startswith('?') and not l.startswith('/'):
                print(" ", l)
except Exception as e:
    print("Error:", e)
