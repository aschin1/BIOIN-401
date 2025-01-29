import requests
import sys
import os
import re
import urllib.parse
from Bio import SeqIO

def get_accessions(query, page_size=500):
    """Retrieves UniProt accessions from UniProtKB based on the given query."""
    link = f"https://rest.uniprot.org/uniprotkb/search?format=json&query={query}&size={page_size}"
    page = 1
    accessions = []
    
    response = requests.get(link)
    if response.status_code != 200:
        print(f"Error: Failed to fetch data from UniProt (Status Code: {response.status_code})")
        print(f"Query used: {link}")  # Debugging
        return []

    data = response.json()
    for entry in data.get("results", []):
        accessions.append(entry["primaryAccession"])

    # Handle pagination
    while "next" in data.get("links", {}):
        next_link = data["links"]["next"]
        response = requests.get(next_link)
        if response.status_code != 200:
            break
        data = response.json()
        for entry in data.get("results", []):
            accessions.append(entry["primaryAccession"])

    return accessions

def extract_sequence(fasta_file):
    """Extracts the sequence from a FASTA file."""
    try:
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                return str(record.seq)
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)
    return None

def get_fasta_metadata(fasta_file):
    """Extracts metadata (Accession, Gene Name) from the FASTA header."""
    try:
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id.split("|")[1] if "|" in record.id else record.id
                gene_name = record.description.split("GN=")[1].split()[0] if "GN=" in record.description else None
                return accession, gene_name
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python accession_retrieval2.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]  
    sequence = extract_sequence(fasta_file)
    accession_number, gene_name = get_fasta_metadata(fasta_file)

    if not sequence:
        print("Error: No sequence found in FASTA file.")
        sys.exit(1)

    print(f"Extracted Accession Number: {accession_number}")
    print(f"Extracted Gene Name: {gene_name}")
    print(f"Extracted Sequence: {sequence}")

    # ✅ Use accession number or gene name for querying UniProt
    if accession_number:
        query = f"accession:{accession_number}"
    elif gene_name:
        query = f"gene:{gene_name}"
    else:
        print("⚠ No valid query parameters found. Exiting.")
        sys.exit(1)

    accessions = get_accessions(query)

    if accessions:
        print("Found accessions:", accessions)
    else:
        print("No accessions found for the given query.")
