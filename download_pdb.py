# Usage: python download_pdb.py <page_from> <page_to>
#
# Script requires `accessions_{page}.txt` files to be present in the same directory as the script.
# You can generate those from get_accessions.py.
#
# You can do thing like `python download_pdb.py 1 10` to download PDB files for accessions from page 1 to 10.
# You can do `python download_pdb.py 1 10 | tee download.log` to log the output and keep track of progress.

import requests
import sys
import os

page_from = sys.argv[1]
page_to = sys.argv[2]
outdir = "pdb_files"

os.system(f"mkdir -p {outdir}")

for page in range(int(page_from), int(page_to)+1):
    with open(f"accessions_{page}.txt", "r") as f:
        accessions = f.readlines()
        print (f"Downloading from page {page}...")
        for accession in accessions:
            accession = accession.strip()
            link = "https://alphafold.ebi.ac.uk/api/prediction/" + accession
            print (link)
            response = requests.get(link)
            pdbUrl = response.json()[0]['pdbUrl']
            os.system(f"wget {pdbUrl} -O ./{outdir}/{page}/{accession}.pdb")