import requests
import sys
import os

def fetch_alphafold_pdb(accession):
    """Fetches PDB structure from AlphaFold for a given UniProt accession."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        print(f"Error: Failed to fetch AlphaFold PDB for {accession} (Status Code: {response.status_code})")
        return None

def save_pdb_file(accession, pdb_data):
    """Saves PDB data to a file in a separate folder."""
    output_dir = "pdb_files"
    os.makedirs(output_dir, exist_ok=True)  # Create directory if it doesn't exist
    
    pdb_filename = os.path.join(output_dir, f"{accession}.pdb")
    with open(pdb_filename, "w") as pdb_file:
        pdb_file.write(pdb_data)
    print(f"PDB file saved as {pdb_filename}")

if __name__ == "__main__":
    accession_number = input("Enter the UniProt accession number: ").strip()
    
    if not accession_number:
        print("Error: No accession number provided.")
        sys.exit(1)
    
    pdb_data = fetch_alphafold_pdb(accession_number)
    
    if pdb_data:
        save_pdb_file(accession_number, pdb_data)
