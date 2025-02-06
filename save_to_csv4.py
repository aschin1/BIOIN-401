import csv
import os

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

