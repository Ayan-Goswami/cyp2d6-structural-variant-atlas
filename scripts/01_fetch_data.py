import pathlib
import requests

files = {
    "structures/2F9Q.pdb": "https://files.rcsb.org/download/2F9Q.pdb",
    "data/raw/P10635.fasta": "https://rest.uniprot.org/uniprotkb/P10635.fasta",
}

for path, url in files.items():
    p = pathlib.Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    p.write_bytes(r.content)
    print(f"Saved {path} ({len(r.content):,} bytes)")