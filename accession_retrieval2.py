import requests
import sys
import os
import re
from Bio import SeqIO

def get_accessions(query, page_size=500):
    """Retrieves UniProt accessions from UniRef based on the given query."""
    link = f"https://rest.uniprot.org/uniref/search?format=list&query={query}&size={page_size}"
    page = 1
    accessions = []
    
    response = requests.get(link)
    if response.status_code != 200:
        print(f"Error: Failed to fetch data from UniProt (Status Code: {response.status_code})")
        return []

    while response.status_code == 200:
        print(f"Getting accession list for page: {page}...")
        lines = response.text.strip().split('\n')

        for line in lines:
            if line:  # Avoid empty lines
                accessions.append(line.strip())

        # Check for next page
        match = re.search(r'<(https://[^>]+)>; rel="next"', response.headers.get('Link', ''))
        next_link = match.group(1) if match else None
        if not next_link:
            break  # Stop if no next page
        
        page += 1
        response = requests.get(next_link)

    return accessions

def extract_sequence(fasta_file):
    """Extracts the sequence from a FASTA file."""
    try:
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                return str(record.seq)  # Return the first sequence found
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)

    return None  # Return None if no sequence is found

def get_accession_number(fasta_file):
    """Extracts the accession number from the FASTA header."""
    try:
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id.split("|")[1] if "|" in record.id else record.id
                return accession
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <fasta_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]  # Extract filename correctly
    sequence = extract_sequence(fasta_file)
    accession_number = get_accession_number(fasta_file)

    if not sequence:
        print("Error: No sequence found in FASTA file.")
        sys.exit(1)

    print(f"Extracted Accession Number: {accession_number}")
    print(f"Extracted Sequence: {sequence}")

    # Construct your UniRef query based on the sequence
    query = f"(sequence:{sequence})"  # Example query, adjust as needed
    accessions = get_accessions(query)

    if accessions:
        print("Found accessions:", accessions)
    else:
        print("No accessions found for the given sequence.")
