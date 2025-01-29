#edited on Jan 27,2025
#extracts sequences from a fasta file and returns them as a list of strings.

from Bio import SeqIO

def extract_sequences(fasta_file):
    """
    Extracts sequences from a FASTA file.

    Args:
        fasta_file: Path to the FASTA file.

    Returns:
        A list of sequences extracted from the FASTA file.
    """
    sequences = []
    with open(fasta_file, "r") as f: 
        for record in SeqIO.parse(f, "fasta"): 
            sequences.append(str(record.seq))
            print(sequences)
    return sequences

sequences = extract_sequences("sample1.txt")
print(f"Number of sequences in file: {len(sequences)}")