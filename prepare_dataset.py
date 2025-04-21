import pandas as pd
import torch
from transformers import AutoTokenizer
from sklearn.model_selection import train_test_split
import csv
import os

# Load dataset
df = pd.read_csv("protein_structures_rand_2000.csv")

# Define label mapping
label_map = {"H": 0, "B": 1, "C": 2}  # 3-class labels

# Load ProtBERT tokenizer
model_name = "Rostlab/prot_bert"
tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)

# === Logging Function ===
def log_filtering_info(seq_id, original_len, filtered_count, reasons, logfile="filtering_log.csv"):
    file_exists = os.path.isfile(logfile)
    with open(logfile, mode="a", newline="") as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(["Sequence ID", "Original Length", "Filtered Tokens", "Reason(s)"])
        writer.writerow([seq_id, original_len, filtered_count, "; ".join(reasons)])

# === Sliding Window Function ===
def sliding_window(sequence, structure, window_size=512, num_segments=10):
    seq_length = len(sequence)
    if seq_length < window_size:
        return [(sequence, structure)]
    step_size = max(1, (seq_length - window_size) // (num_segments - 1))
    return [
        (sequence[i : i + window_size], structure[i : i + window_size]) 
        for i in range(0, seq_length - window_size + 1, step_size)
    ][:num_segments]

# === Tokenization + Logging ===
def tokenize_and_align_labels(sequence, structure, seq_id=None, max_len=512):
    formatted_sequence = " ".join(sequence)
    structure_labels = [label_map.get(res, -100) for res in structure]
    tokenized_inputs = tokenizer(formatted_sequence, padding="max_length", truncation=True, max_length=max_len)
    attention_mask = tokenized_inputs["attention_mask"]
    aligned_labels = [-100] * len(attention_mask)

    seq_index = 0
    for j, mask in enumerate(attention_mask):
        if mask == 1 and seq_index < len(structure_labels):
            aligned_labels[j] = structure_labels[seq_index]
            seq_index += 1

    # === Filtering Log Logic ===
    filtered = aligned_labels.count(-100)
    total = len(structure)
    reasons = []
    if total != sum(mask == 1 for mask in attention_mask):
        reasons.append("label/token mismatch")
    if filtered >= len(aligned_labels):
        reasons.append("all tokens filtered")
    if -100 in structure_labels:
        reasons.append("unknown structure label")

    if seq_id is not None:
        log_filtering_info(seq_id, total, filtered, reasons)

    return dict(tokenized_inputs), torch.tensor(aligned_labels)

# === Process dataset with logging ===
tokenized_sequences = []
label_sequences = []

for idx, row in df.iterrows():
    segments = sliding_window(row["Primary Sequence"], row["Secondary Structure"])
    for segment_id, (seq_chunk, label_chunk) in enumerate(segments):
        full_id = f"Row{idx}_Seg{segment_id}"
        chunk_inputs, chunk_labels = tokenize_and_align_labels(seq_chunk, label_chunk, seq_id=full_id)
        tokenized_sequences.append(chunk_inputs)
        label_sequences.append(chunk_labels)

# === Split and Save ===
train_inputs, test_inputs, train_labels, test_labels = train_test_split(
    tokenized_sequences, label_sequences, test_size=0.2, random_state=42
)

torch.save((train_inputs, train_labels), "train_data.pt")
torch.save((test_inputs, test_labels), "test_data.pt")

print("✅ Data prepared for ProtBERT fine-tuning with Sliding Window!")
print("📝 Filtering log saved to filtering_log.csv")
