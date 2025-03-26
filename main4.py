import sys
import os
import csv
from accession_retrieval4 import get_accessions, extract_sequence, get_fasta_metadata, get_accession_by_sequence
from get_alphafold_data4 import fetch_alphafold_pdb, save_pdb_file
from get_secondary_struct4 import run_phipsi, run_define2, parse_define2_output, display_horizontal_output

def check_existing_sequence(sequence, dataset_file="protein_structures_3500.csv"):
    """
    Checks if a given sequence exists in the dataset and returns the corresponding secondary structure.
    """
    if not os.path.isfile(dataset_file):
        return None
    
    with open(dataset_file, mode="r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["Primary Sequence"] == sequence:
                print("✅ Sequence found in dataset! Returning existing secondary structure.")
                return row["Secondary Structure"]
    return None

def save_to_csv(accession, primary_sequence, secondary_structure, output_file="protein_structures_3500.csv"):
    """
    Saves the accession number, primary sequence, and secondary structure labels to a CSV file.
    """
    file_exists = os.path.isfile(output_file)
    with open(output_file, mode="a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        if not file_exists:
            writer.writerow(["Accession Number", "Primary Sequence", "Secondary Structure"])
        writer.writerow([accession, primary_sequence, secondary_structure])
    print(f"Data saved to {output_file}")

def main():
    if len(sys.argv) == 2:
        fasta_file = sys.argv[1]
        if not os.path.isfile(fasta_file):
            print(f"Error: File '{fasta_file}' not found.")
            sys.exit(1)
        sequence = extract_sequence(fasta_file)
        accession_number, gene_name = get_fasta_metadata(fasta_file)
    else:
        option = input("Choose an option:\n1. Input a FASTA file\n2. Input a protein sequence directly\nEnter 1 or 2: ")
        if option == "1":
            fasta_file = input("Enter the path to your FASTA file: ")
            if not os.path.isfile(fasta_file):
                print(f"Error: File '{fasta_file}' not found.")
                sys.exit(1)
            sequence = extract_sequence(fasta_file)
            accession_number, gene_name = get_fasta_metadata(fasta_file)
        elif option == "2":
            sequence = input("Enter your protein sequence: ").strip()
            if not sequence:
                print("Error: No sequence provided.")
                sys.exit(1)
            accession_number, gene_name = None, None
            accessions = get_accession_by_sequence(sequence)
            if accessions:
                accession_number = accessions[0]
        else:
            print("Invalid option. Exiting.")
            sys.exit(1)
    
    # Check if the secondary structure is already known
    known_secondary_structure = check_existing_sequence(sequence)
    if known_secondary_structure:
        print("Secondary Structure:", known_secondary_structure)
        return
    
    if not accession_number:
        print("Error: No accession number found.")
        sys.exit(1)
    
    # Fetch the AlphaFold PDB file
    pdb_data = fetch_alphafold_pdb(accession_number)
    if not pdb_data:
        print("Error: Could not fetch PDB file.")
        sys.exit(1)
    
    # Save the PDB file
    pdb_filename = os.path.join("pdb_files", f"{accession_number}.pdb")
    save_pdb_file(accession_number, pdb_data)
    
    # Run phipsi and define2 to get secondary structure
    phipsi_output = run_phipsi(pdb_filename)
    if not phipsi_output:
        print("Error: Phipsi execution failed.")
        sys.exit(1)
    
    define2_output = run_define2(phipsi_output)
    if not define2_output:
        print("Error: Define2 execution failed.")
        sys.exit(1)
    
    residues, sec_structure = parse_define2_output(define2_output)
    if not residues or not sec_structure:
        print("Error: No secondary structure data extracted.")
        sys.exit(1)
    
    # Display secondary structure output
    display_horizontal_output(residues, sec_structure)
    
    # Save results to CSV
    save_to_csv(accession_number, "".join(residues), "".join(sec_structure))

if __name__ == "__main__":
    main()
