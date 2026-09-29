import urllib.request
import re

url = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12063144/'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as resp:
        html = resp.read().decode('utf-8')
        links = re.findall(r'href=[\'\"]([^\'\"]*mmc[^\'\"]*)[\'\"]', html)
        print('Found mmc links:', set(links))
except Exception as e:
    print('Error:', e)
