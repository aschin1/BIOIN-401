import requests
from Bio.PDB import PDBParser, DSSP
import io
import subprocess
import sys
import tempfile

def fetch_alphafold_pdb(accession):
    """Fetches PDB file from AlphaFold for a given accession number."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
    response = requests.get(url)
    print(response.text)
    return response.text if response.status_code == 200 else None

def extract_secondary_structure(pdb_text):
    """Extracts secondary structure sequence labels from a PDB file using DSSP."""
    pdb_parser = PDBParser(QUIET=True)
    
    # 🔹 Create a temporary PDB file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb:
        temp_pdb.write(pdb_text.encode('utf-8'))
        temp_pdb_path = temp_pdb.name  # Store file path

    try:
        structure = pdb_parser.get_structure("AF", temp_pdb_path)
        model = structure[0]
        
        # 🔹 Run DSSP on the actual file (not StringIO)
        dssp = DSSP(model, temp_pdb_path)

        # Extract secondary structure sequence
        secondary_structure = "".join([dssp[key][2] for key in dssp.keys()])
        
    finally:
        # 🔹 Remove the temporary file after processing
        os.remove(temp_pdb_path)
    
    return secondary_structure

def get_accession_from_fasta(fasta_file):
    """Runs accession_retrieval2.py to extract accession number from a FASTA file."""
    result = subprocess.run(["python", "accession_retrieval2.py", fasta_file], capture_output=True, text=True)
    
    for line in result.stdout.split("\n"):
        if "Extracted Accession Number:" in line:
            return line.split(":")[1].strip()
    
    return None

def process_fasta(fasta_file):
    accession = get_accession_from_fasta(fasta_file)
    
    if not accession:
        print("Error: No valid accession number found.")
        return None
    
    print(f"Using Accession: {accession} for AlphaFold PDB retrieval.")
    pdb_text = fetch_alphafold_pdb(accession)
    
    if pdb_text:
        return extract_secondary_structure(pdb_text)
    
    print(f"Failed to retrieve PDB data for accession: {accession}")
    return None

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Extract secondary structure labels from AlphaFold PDB files using an accession number retrieved from a FASTA file.")
    parser.add_argument("fasta_file", type=str, help="Path to the FASTA file to retrieve accession number and fetch PDB structure.")
    
    args = parser.parse_args()
    structure_labels = process_fasta(args.fasta_file)
    
    if structure_labels:
        print(f"Secondary Structure: {structure_labels}")
    else:
        print("Failed to retrieve or process PDB data.")
