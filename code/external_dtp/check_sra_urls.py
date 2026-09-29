import urllib.request
import re

url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=43550465&retmode=xml'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as resp:
        xml = resp.read().decode('utf-8')
        urls = re.findall(r'url=[\'"]([^\'"]+)[\'"]', xml)
        print('URLs in SRA XML:')
        for u in set(urls):
            print(' ', u)
        sizes = re.findall(r'size=[\'"]([^\'"]+)[\'"]', xml)
        print('Sizes:', set(sizes))
except Exception as e:
    print('Error:', e)
