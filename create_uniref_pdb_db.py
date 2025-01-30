#! /usr/bin/env python3

import argparse
import requests
import subprocess
import pandas as pd
from proteinbert.shared_utils.util import get_parser_file_type
from proteinbert.uniref_dataset import UnirefToSqliteParser

def get_accession_from_fasta(fasta_file):
    """Runs accession_retrieval2.py to extract accession number from a FASTA file."""
    result = subprocess.run(["python", "accession_retrieval2.py", fasta_file], capture_output=True, text=True)
   # print(result)
    for line in result.stdout.split("\n"):
        if "Extracted Accession Number:" in line:
            return line.split(":")[1].strip()
    
    return None

def fetch_alphafold_pdb(accession):
    """Fetches PDB structure from AlphaFold for a given UniProt accession."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
    response = requests.get(url)
    return response.text if response.status_code == 200 else None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create an SQLite DB with AlphaFold PDB structures from a raw UniRef file.')
    parser.add_argument('--fasta-file', dest='fasta_file', metavar='/path/to/sequence.fasta', 
                        type=get_parser_file_type(parser, must_exist=True), required=True, help='Path to the FASTA file.')
    parser.add_argument('--output-sqlite-file', dest='output_sqlite_file', metavar='/path/to/uniref.db', 
                        type=get_parser_file_type(parser), required=True, help='Path to save the output SQLite file.')
    parser.add_argument('--log-progress-every', dest='log_progress_every', metavar='1000', type=int, default=1000, 
                        help='Log progress in increments of this given number (default 1000).')
    parser.add_argument('--chunk-size', dest='chunk_size', metavar='100000', type=int, default=100000, 
                        help='Number of protein records per chunk written into the created DB.')
    parser.add_argument('--silent', dest='silent', action='store_true', help='Run in silent mode.')
    args = parser.parse_args()
    
    # Retrieve accession number from FASTA file
    accession = get_accession_from_fasta(args.fasta_file)
    print(accession)
    if not accession:
        print("Error: No valid accession number found.")
        exit(1)
    
    # Fetch PDB structure using the retrieved accession number
    pdb_structure = fetch_alphafold_pdb(accession)
    
    # Initialize SQLite parser
    parser = UnirefToSqliteParser(None, None, args.output_sqlite_file, 
                                  verbose=not args.silent, log_progress_every=args.log_progress_every, 
                                  chunk_size=args.chunk_size)
    
    # Store extracted PDB data in SQLite
    record = {"accession": accession, "pdb_structure": pdb_structure}
    parser.store(record)
    
    print(f"Database saved at {args.output_sqlite_file}")
