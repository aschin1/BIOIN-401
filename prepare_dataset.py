import pandas as pd
import torch
from transformers import AutoTokenizer
from sklearn.model_selection import train_test_split

# Load dataset
df = pd.read_csv("protein_structures_50.csv")

# Define label mapping
label_map = {"H": 0, "B": 1, "C": 2}  # 3-class labels

# Load ProtBERT tokenizer
model_name = "Rostlab/prot_bert"
tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)

def sliding_window(sequence, structure, window_size=512, num_segments=10):
    """
    Generate overlapping segments of a protein sequence and its secondary structure labels.

    :param sequence: Protein sequence (string).
    :param structure: Corresponding secondary structure (string).
    :param window_size: Length of each segment (default: 512).
    :param num_segments: Number of segments to generate (default: 10).
    :return: List of (sequence segments, label segments).
    """
    seq_length = len(sequence)
    
    if seq_length < window_size:
        return [(sequence, structure)]  # If sequence is too short, return as a single segment

    # Compute step size for even overlap
    step_size = max(1, (seq_length - window_size) // (num_segments - 1))

    segments = [
        (sequence[i : i + window_size], structure[i : i + window_size]) 
        for i in range(0, seq_length - window_size + 1, step_size)
    ]

    return segments[:num_segments]  # Ensure exactly `num_segments` windows

def tokenize_and_align_labels(sequence, structure, max_len=512):
    """
    Tokenizes protein sequence and aligns secondary structure labels with tokenized output.

    :param sequence: Protein sequence segment (string).
    :param structure: Corresponding secondary structure (string).
    :param max_len: Max sequence length (default: 512).
    :return: Tokenized inputs and aligned labels.
    """
    formatted_sequence = " ".join(sequence)  # Space-separated residues for tokenization
    structure_labels = [label_map.get(res, -100) for res in structure]  # Convert to numerical labels

    tokenized_inputs = tokenizer(formatted_sequence, padding="max_length", truncation=True, max_length=max_len)

    # Align labels with tokenized input
    attention_mask = tokenized_inputs["attention_mask"]
    aligned_labels = [-100] * len(attention_mask)  # Initialize ignored tokens

    seq_index = 0
    for j, mask in enumerate(attention_mask):
        if mask == 1 and seq_index < len(structure_labels):
            aligned_labels[j] = structure_labels[seq_index]
            seq_index += 1

    return dict(tokenized_inputs), torch.tensor(aligned_labels)

# Process dataset with sliding window
tokenized_sequences = []
label_sequences = []

for _, row in df.iterrows():
    segments = sliding_window(row["Primary Sequence"], row["Secondary Structure"])  # Apply sliding window
    for seq_chunk, label_chunk in segments:
        chunk_inputs, chunk_labels = tokenize_and_align_labels(seq_chunk, label_chunk)
        tokenized_sequences.append(chunk_inputs)
        label_sequences.append(chunk_labels)

# Split dataset into train and test
train_inputs, test_inputs, train_labels, test_labels = train_test_split(
    tokenized_sequences, label_sequences, test_size=0.2, random_state=42
)

# Save preprocessed data
torch.save((train_inputs, train_labels), "train_data.pt")
torch.save((test_inputs, test_labels), "test_data.pt")

print("✅ Data prepared for ProtBERT fine-tuning with Sliding Window!")
