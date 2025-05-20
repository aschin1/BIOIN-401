from Bio import SeqIO
import os
import subprocess
import tempfile

def run_iupred2(seq, mode="long", threshold=0.5):
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f:
        f.write(">query\n" + seq)
        fasta_path = f.name

    output_path = fasta_path + ".out"
    try:
        subprocess.run(["python", "iupred2a/iupred2a/iupred2a.py", fasta_path, mode],
    stdout=open(output_path, "w"),
    check=True
)


        return parse_iupred_output(output_path, threshold)
    finally:
        os.remove(fasta_path)
        if os.path.exists(output_path):
            os.remove(output_path)

def parse_iupred_output(path, threshold=0.5):
    disordered = 0
    total = 0
    regions = []
    current = []

    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            pos, aa, score = parts
            score = float(score)
            total += 1
            if score > threshold:
                disordered += 1
                current.append(int(pos))
            else:
                if current:
                    regions.append((current[0], current[-1]))
                    current = []

    if current:
        regions.append((current[0], current[-1]))

    fraction = disordered / total if total else 0
    classification = "disordered" if fraction > 0.5 else "structured"

    return {
        "classification": classification,
        "disorder_fraction": round(fraction, 3),
        "regions": regions
    }

# Usage example
if __name__ == "__main__":
    fasta_file = input("Enter the path to your FASTA file: ").strip()
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq = str(record.seq)
            result = run_iupred2(seq)
            print(f"Classification: {result['classification']}")
            print(f"Disorder Fraction: {result['disorder_fraction']}")
            print("Disordered Regions:", result['regions'])
