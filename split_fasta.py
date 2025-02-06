import os
from Bio import SeqIO

def split_fasta(input_fasta, output_folder="primary_sequences"):
    """
    Splits a multi-entry FASTA file into individual FASTA files.
    
    :param input_fasta: Path to the input FASTA file.
    :param output_folder: Directory where individual FASTA files will be stored.
    """
    # Ensure the output directory exists
    os.makedirs(output_folder, exist_ok=True)

    # Read and process each sequence in the FASTA file
    with open(input_fasta, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract the identifier (use the first part of the header)
            fasta_id = record.id.split()[0]

            # Define output file path
            output_fasta_path = os.path.join(output_folder, f"{fasta_id}.fasta")

            # Write the sequence to an individual FASTA file
            with open(output_fasta_path, "w") as output_file:
                SeqIO.write(record, output_file, "fasta")

            print(f"Saved: {output_fasta_path}")

    print("\n✅ Splitting completed! All sequences are saved in the 'primary_sequences' folder.")

# Example usage:
if __name__ == "__main__":
    input_fasta_file = input("Enter the path to the multi-entry FASTA file: ").strip()
    split_fasta(input_fasta_file)
