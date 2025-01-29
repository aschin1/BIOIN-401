# Usage: python get_accessions.py
#
# Note that because of how UniProt's API pagination works, we can't specify number of pages to download nor which
# page to start from. You'll have to make sure that script runs to completion to get all the accessions.

import requests
import sys
import os


link = "https://rest.uniprot.org/uniref/search?format=list&query=%28%28identity%3A0.5%29%29&size=500"
page = 1
response = requests.get(link)

while response.headers['Link']:
    print ("Getting accession list for page: " + str(page) + "...")
    lines = response.text.split('\n')
    accessions = []

    for link in lines:
        accession = link[9:]
        os.system(f"echo {accession} >> accessions_{page}.txt")
        
    next_link = response.headers['Link'].split('<')[1].split('>; rel="next"')[0]
    page += 1
    response = requests.get(next_link)