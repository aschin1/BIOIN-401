import os
import random
import shutil

# Directory containing the 572,000 FASTA files
SOURCE_DIR = "transmembrane_protein_sequences"  # Change if your FASTA files are elsewhere

# Desired output sample sizes and corresponding output folders
SAMPLE_SIZES = [1000, 2000]

def create_randomized_sets(source_dir, sample_sizes):
    # Get list of all FASTA files in source directory
    fasta_files = [f for f in os.listdir(source_dir) if f.endswith(".fasta")]
    print(f"Found {len(fasta_files)} FASTA files.")

    # Shuffle once to maintain consistent randomness across all samples
    random.shuffle(fasta_files)

    for size in sample_sizes:
        sample = fasta_files[:size]
        output_dir = f"transmembrane_subset_n{size}"
        os.makedirs(output_dir, exist_ok=True)

        for file_name in sample:
            source_path = os.path.join(source_dir, file_name)
            dest_path = os.path.join(output_dir, file_name)
            shutil.copy2(source_path, dest_path)

        print(f"✅ Created folder '{output_dir}' with {size} randomized fasta files.")

if __name__ == "__main__":
    create_randomized_sets(SOURCE_DIR, SAMPLE_SIZES)
