import sys
import os
import csv
from accession_retrieval4 import extract_sequence
from get_secondary_struct_dssp import get_dssp_secondary_structure, display_horizontal_output as display_dssp
from get_secondary_struct4 import run_phipsi, run_define2, parse_define2_output, display_horizontal_output as display_vadar

def is_entry_in_dataset(accession, dataset_file="pdb_new.csv"):
    if not os.path.isfile(dataset_file):
        return False
    with open(dataset_file, mode="r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row["Accession Number"] == accession:
                return True
    return False

def save_to_csv(accession, primary_sequence, vadar_structure, dssp_structure, output_file="pdb_new.csv"):
    file_exists = os.path.isfile(output_file)
    with open(output_file, mode="a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        if not file_exists:
            writer.writerow(["Accession Number", "Primary Sequence", "VADAR Structure", "DSSP Structure"])
        writer.writerow([accession, primary_sequence, vadar_structure, dssp_structure])
    print(f"✅ Data saved to {output_file}")

def compute_similarity(vadar_structure, dssp_structure):
    if len(vadar_structure) != len(dssp_structure):
        print("Error: Length mismatch between VADAR and DSSP structures.")
        return None
    similarity = sum(1 for v, d in zip(vadar_structure, dssp_structure) if v == d) / len(vadar_structure) * 100
    print(f"Similarity: {similarity:.2f}%")
    return similarity

def main():
    if len(sys.argv) == 2:
        pdb_file = sys.argv[1]
        if not os.path.isfile(pdb_file):
            print(f"Error: File '{pdb_file}' not found.")
            sys.exit(1)
        accession_number = os.path.splitext(os.path.basename(pdb_file))[0]
    else:
        pdb_file = input("Enter the path to your PDB file: ")
        if not os.path.isfile(pdb_file):
            print(f"Error: File '{pdb_file}' not found.")
            sys.exit(1)
        accession_number = os.path.splitext(os.path.basename(pdb_file))[0]

    if is_entry_in_dataset(accession_number):
        print(f"✅ Entry for accession {accession_number} already exists in the dataset. Skipping.")
        return

    phipsi_output = run_phipsi(pdb_file)
    if not phipsi_output:
        print("Error: Phipsi execution failed.")
        sys.exit(1)

    define2_output = run_define2(phipsi_output)
    if not define2_output:
        print("Error: Define2 execution failed.")
        sys.exit(1)

    vadar_residues, vadar_structure = parse_define2_output(define2_output)
    if not vadar_residues or not vadar_structure:
        print("Error: No secondary structure data extracted from VADAR.")
        sys.exit(1)

    dssp_residues, dssp_structure = get_dssp_secondary_structure(pdb_file)
    if not dssp_residues or not dssp_structure:
        print("Error: DSSP failed to extract secondary structure.")
        sys.exit(1)

    display_vadar(vadar_residues, vadar_structure)
    display_dssp(dssp_residues, dssp_structure)

    if len(vadar_residues) != len(dssp_residues):
        print("Error: Length mismatch between VADAR and DSSP structures.")
        sys.exit(1)

    similarity = compute_similarity(vadar_structure, dssp_structure)
    if similarity is not None:
        print(f"Similarity between VADAR and DSSP structures: {similarity:.2f}%")

    save_to_csv(accession_number, ''.join(dssp_residues), ''.join(vadar_structure), ''.join(dssp_structure))

if __name__ == "__main__":
    main()
