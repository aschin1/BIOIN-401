import requests
import sys
import time
import urllib.parse
from Bio import SeqIO

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

def get_accession_by_sequence(sequence):
    """
    Submits a protein sequence to EBI's BLAST API to find the closest UniProtKB accession.
    This function polls the API until a result is available.
    """
    blast_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
    check_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/"
    result_url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/"
    
    # Parameters for BLAST search against UniProtKB
    params = {
        "email": "your-email@example.com",  # REQUIRED: Change this to your email
        "program": "blastp",
        "database": "uniprotkb",
        "sequence": sequence
    }
    
    # Submit the BLAST job
    response = requests.post(blast_url, data=params)
    
    if response.status_code != 200:
        print(f"Error: Failed to submit BLAST request (Status Code: {response.status_code})")
        return None

    job_id = response.text.strip()
    print(f"BLAST Job ID: {job_id}")

    # Polling to check the status of the job
    for _ in range(15):  # Wait up to 45 seconds
        status_response = requests.get(check_url + job_id)
        if status_response.text.strip() == "FINISHED":
            break
        print("Waiting for BLAST results...")
        time.sleep(3)

    # Fetch the result if job is finished
    result_response = requests.get(result_url + job_id + "/xml")
    
    if result_response.status_code != 200:
        print(f"Error: Failed to retrieve BLAST results (Status Code: {result_response.status_code})")
        return None

    # Extract the best UniProt accession from the BLAST XML output
    import xml.etree.ElementTree as ET
    root = ET.fromstring(result_response.text)
    
    # Look for hits in the BLAST XML output
    for hit in root.findall(".//Hit"):
        accession = hit.find("Hit_accession").text
        if accession:
            return accession  # Return the first (best) hit

    print("Error: No matching accession found for the given sequence.")
    return None

if __name__ == "__main__":
    option = input("Choose an option:\n1. Input a FASTA file\n2. Input a protein sequence directly\nEnter 1 or 2: ")

    if option == "1":
        fasta_file = input("Enter the path to your FASTA file: ")
        sequence = extract_sequence(fasta_file)
        accession_number, gene_name = get_fasta_metadata(fasta_file)
    elif option == "2":
        sequence = input("Enter your protein sequence: ").strip()
        accession_number, gene_name = None, None
        if not sequence:
            print("Error: No sequence provided.")
            sys.exit(1)
        accessions = get_accession_by_sequence(sequence)
    else:
        print("Invalid option. Exiting.")
        sys.exit(1)

    if option == "1":
        query = f"accession:{accession_number}" if accession_number else f"sequence:{sequence}"
        accessions = get_accessions(query)

    if accessions:
        print("Found accessions:", accessions)
    else:
        print("No accessions found for the given query.")
