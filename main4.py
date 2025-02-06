import sys
import os
from accession_retrieval4 import get_accessions, extract_sequence, get_fasta_metadata, get_accession_by_sequence
from get_alphafold_data4 import fetch_alphafold_pdb, save_pdb_file
from get_secondary_struct4 import run_stride, parse_stride_output, display_horizontal_output
import csv

def save_to_csv(accession, primary_sequence, secondary_structure, output_file="protein_structures.csv"):
    """
    Saves the accession number, primary sequence, and secondary structure labels to a CSV file.

    :param accession: The UniProt accession number.
    :param primary_sequence: The protein sequence.
    :param secondary_structure: The secondary structure labels.
    :param output_file: The name of the CSV file (default: protein_structures.csv).
    """
    file_exists = os.path.isfile(output_file)

    with open(output_file, mode="a", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # Write headers only if the file is newly created
        if not file_exists:
            writer.writerow(["Accession Number", "Primary Sequence", "Secondary Structure"])

        # Write the data
        writer.writerow([accession, primary_sequence, secondary_structure])

    print(f"Data saved to {output_file}")

def main():
    option = input("Choose an option:\n1. Input a FASTA file\n2. Input a protein sequence directly\nEnter 1 or 2: ")

    if option == "1":
        fasta_file = input("Enter the path to your FASTA file: ")
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
            accession_number = accessions[0]  # Use the first accession found
    else:
        print("Invalid option. Exiting.")
        sys.exit(1)

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

    # Run STRIDE to get secondary structure
    stride_output = run_stride(pdb_filename)
    if not stride_output:
        print("Error: STRIDE execution failed.")
        sys.exit(1)

    residues, sec_structure = parse_stride_output(stride_output)
    if not residues or not sec_structure:
        print("Error: No secondary structure data extracted.")
        sys.exit(1)

    # Display secondary structure output
    display_horizontal_output(residues, sec_structure)

    # Save results to CSV
    save_to_csv(accession_number, "".join(residues), "".join(sec_structure))

if __name__ == "__main__":
    main()
