import sys
import os
from Bio import SeqIO
from classify_disorder import run_iupred2, parse_iupred_output
from classify_membrane import run_deeptmhmm  

def read_sequence_input(source):
    if os.path.isfile(source):
        # From FASTA
        for record in SeqIO.parse(source, "fasta"):
            return str(record.seq)
    else:
        return source.strip().upper()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python main_classifier.py <sequence_or_fasta_file>")
        sys.exit(1)

    seq_input = read_sequence_input(sys.argv[1])

    # Disorder classification
    disorder_result = run_iupred2(seq_input)

    print(f"\n🔍 Disorder Classification: {disorder_result['classification']}")
    print(f"  Disorder fraction: {disorder_result['disorder_fraction']}")
    print(f"  Disordered regions: {disorder_result['regions']}")

    # Membrane classification
    membrane_class, tm_regions = run_deeptmhmm(seq_input)
    print(f"\n🔍 Membrane Classification: {membrane_class}")
    if tm_regions:
        print(f"  TM Helices: {tm_regions}")

    # 🧬 Final combined classification
    print("\n=== Final Protein Classification ===")
    if membrane_class == "membrane":
        final_class = "membrane"
    elif disorder_result["classification"] == "disordered":
        final_class = "disordered"
    else:
        final_class = "structured"

    print(f"🧬 Protein classified as: {final_class}")

    

