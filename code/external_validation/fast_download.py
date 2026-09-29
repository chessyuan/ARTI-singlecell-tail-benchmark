import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4150nnn/GSM4150377/suppl/GSM4150377_sciPlex2_A549_Transcription_Modulators_UMI.count.matrix.gz"
dest = "data_raw/external_acute/sciplex2/GSM4150377_sciPlex2_A549_Transcription_Modulators_UMI.count.matrix.gz"
part_dir = dest + ".parts"
os.makedirs(part_dir, exist_ok=True)

proxy_handler = urllib.request.ProxyHandler({'http': 'http://127.0.0.1:7897', 'https': 'http://127.0.0.1:7897'})
opener = urllib.request.build_opener(proxy_handler)

req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': 'Mozilla/5.0'})
with opener.open(req, timeout=15) as resp:
    total_size = int(resp.headers.get('Content-Length', 0))

print(f"Total file size: {total_size} bytes ({total_size/(1024*1024):.2f} MB)")

num_chunks = 16
chunk_size = (total_size + num_chunks - 1) // num_chunks

def download_chunk(idx):
    start = idx * chunk_size
    end = min(total_size - 1, (idx + 1) * chunk_size - 1)
    if start > end:
        return
    part_file = os.path.join(part_dir, f"part_{idx:02d}")
    if os.path.exists(part_file) and os.path.getsize(part_file) == (end - start + 1):
        # already done
        return
    
    # download with retry
    for attempt in range(5):
        try:
            r = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0',
                'Range': f'bytes={start}-{end}'
            })
            with opener.open(r, timeout=30) as res:
                data = res.read()
                if len(data) == (end - start + 1):
                    with open(part_file, 'wb') as pf:
                        pf.write(data)
                    print(f"Chunk {idx:02d} completed ({len(data)} bytes).")
                    return
                else:
                    print(f"Chunk {idx:02d} incomplete ({len(data)} / {end - start + 1}), retrying...")
        except Exception as e:
            print(f"Chunk {idx:02d} attempt {attempt+1} failed: {e}")
            time.sleep(2)
    raise RuntimeError(f"Chunk {idx:02d} failed after 5 attempts.")

print(f"Starting {num_chunks} parallel chunk downloads...")
t0 = time.time()
with ThreadPoolExecutor(max_workers=16) as pool:
    list(pool.map(download_chunk, range(num_chunks)))
t1 = time.time()
print(f"All chunks downloaded in {t1 - t0:.2f}s! Merging...")

with open(dest, 'wb') as outfile:
    for idx in range(num_chunks):
        part_file = os.path.join(part_dir, f"part_{idx:02d}")
        with open(part_file, 'rb') as pf:
            outfile.write(pf.read())
        os.remove(part_file)
os.rmdir(part_dir)
print(f"Successfully assembled {dest}: {os.path.getsize(dest)} bytes!")
