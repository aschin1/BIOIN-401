import requests
import sys
import time
import urllib.parse
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def get_accessions(query, page_size=500):
    """Retrieves UniProt accessions from UniProtKB based on the given query."""
    query = urllib.parse.quote(query)  # URL encode query
    link = f"https://rest.uniprot.org/uniprotkb/search?format=json&query={query}&size={page_size}"
    print("Requesting URL:", link)  # Debugging print
    accessions = []

    response = requests.get(link)
    print("Response Status Code:", response.status_code)  # Debugging print
    if response.status_code != 200:
        print(f"Error: Failed to fetch data from UniProt (Status Code: {response.status_code})")
        return []

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

def get_accession_by_sequence(seq, identity_threshold=95.0):
    print("Running BLAST against UniProtKB (SwissProt)...")
    result_handle = NCBIWWW.qblast("blastp", "swissprot", seq)

    blast_record = NCBIXML.read(result_handle)

    if not blast_record.alignments:
        print("❌ No matches found.")
        return None, None

    first_alignment = blast_record.alignments[0]
    best_hsp = first_alignment.hsps[0]

    accession = first_alignment.accession
    identity = (best_hsp.identities / best_hsp.align_length) * 100

    print(f"🔍 Best hit: {accession} ({identity:.2f}% identity)")

    if identity >= identity_threshold:
        print(f"✅ Identity above threshold ({identity_threshold}%) — returning accession.")
        return accession, identity
    else:
        print(f"❌ Identity {identity:.2f}% is below threshold ({identity_threshold}%).")
        return None, None

if __name__ == "__main__":
    option = input("Choose an option:\n1. Input a FASTA file\n2. Input a protein sequence directly\nEnter 1 or 2: ")

    if option == "1":
        fasta_file = input("Enter the path to your FASTA file: ")
        sequence = extract_sequence(fasta_file)
        accession_number, gene_name = get_fasta_metadata(fasta_file)
        query = f"accession:{accession_number}" if accession_number else f"sequence:{sequence}"
        accessions = get_accessions(query)

    elif option == "2":
        sequence = input("Enter your protein sequence: ").strip()
        accession_number, gene_name = None, None
        if not sequence:
            print("Error: No sequence provided.")
            sys.exit(1)
        accession_number, identity = get_accession_by_sequence(sequence, identity_threshold=95.0)
        accessions = [accession_number] if accession_number else []

    else:
        print("Invalid option. Exiting.")
        sys.exit(1)

    if accessions:
        print("✅ Found accessions:", accessions)
    else:
        print("❌ No accessions found for the given query.")