import os
import tempfile
import subprocess
from Bio import SeqIO

def run_deeptmhmm(seq):
    """
    Runs DeepTMHMM via Biolib CLI and checks for transmembrane prediction.
    Returns:
        - classification: 'membrane' or 'non-membrane'
        - tm_regions: list of residue ranges with TM annotation
    """
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f:
        f.write(">query\n" + seq)
        fasta_path = f.name

    try:
        # Run DeepTMHMM through Biolib
        subprocess.run(["biolib", "run", "DTU/DeepTMHMM", "--fasta", fasta_path], check=True)

        # Output file location (by default)
        result_path = "biolib_results/deeptmhmm_results.md"
        if not os.path.exists(result_path):
            print("❌ DeepTMHMM result file not found.")
            return "error", []

        tm_regions = []
        with open(result_path, "r") as f:
            for line in f:
                if "TMhelix" in line or "TM" in line:
                    tokens = line.strip().split()
                    # Try to extract positions if available
                    for i, t in enumerate(tokens):
                        if t.startswith("TM") and i + 2 < len(tokens):
                            try:
                                start = int(tokens[i+1])
                                end = int(tokens[i+2])
                                tm_regions.append((start, end))
                            except ValueError:
                                continue

        classification = "membrane" if tm_regions else "non-membrane"
        return classification, tm_regions

    finally:
        os.remove(fasta_path)
if __name__ == "__main__":
    # Example usage
    fasta_file = input("Enter the path to your FASTA file: ").strip()
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq = str(record.seq)
            classification, tm_regions = run_deeptmhmm(seq)
            print(f"Classification: {classification}")
            if tm_regions:
                print("TM Regions:", tm_regions)
            else:
                print("No TM regions found.")