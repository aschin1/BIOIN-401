import torch
from transformers import AutoTokenizer
from Bio import pairwise2
from tqdm import tqdm

def decode_sequence(input_ids, tokenizer):
    """Convert token IDs to amino acid letters, skipping special tokens."""
    tokens = tokenizer.convert_ids_to_tokens(input_ids)
    amino_acids = [tok for tok in tokens if tok not in tokenizer.all_special_tokens]
    return "".join(amino_acids)

def full_sequence_identity(seq1, seq2):
    """Compute % identity using global alignment."""
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    if not alignments:
        return 0.0
    aligned_seq1, aligned_seq2, _, _, _ = alignments[0]
    matches = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2))
    return (matches / max(len(aligned_seq1), 1)) * 100

# Load tokenized datasets
train_inputs, _ = torch.load("train.pt")
val_inputs, _ = torch.load("validation.pt")

# Load tokenizer
tokenizer = AutoTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)

# Decode sequences
train_seqs = [decode_sequence(sample["input_ids"], tokenizer) for sample in train_inputs]
val_seqs = [decode_sequence(sample["input_ids"], tokenizer) for sample in val_inputs]

# ✅ Detect and report exact overlaps
print("\n🔍 Checking for exact overlapping sequences...\n")
val_seq_set = set(val_seqs)
exact_overlaps = [seq for seq in train_seqs if seq in val_seq_set]

if exact_overlaps:
    print(f"❌ Found {len(exact_overlaps)} exact overlapping sequences between training and validation sets!")
else:
    print("✅ No exact sequence overlaps found.")

# ✅ Compute average % identity (excluding exact matches)
print("\n🔬 Comparing training vs. validation sequence identity...\n")

total_identity = 0.0
pair_count = 0

for train_seq in tqdm(train_seqs, desc="Training Sequences"):
    for val_seq in val_seqs:
        if train_seq == val_seq:
            continue  # Skip exact matches
        identity = full_sequence_identity(train_seq, val_seq)
        total_identity += identity
        pair_count += 1

# Final reporting
if pair_count > 0:
    avg_identity = total_identity / pair_count
    print(f"\n✅ Average % sequence identity between train and validation (non-identical): {avg_identity:.2f}%")
    print(f"🔢 Compared {pair_count} unique sequence pairs.")
else:
    print("⚠️ No unique train/val pairs found for comparison.")
