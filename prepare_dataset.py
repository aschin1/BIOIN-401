import pandas as pd
import torch
from transformers import AutoTokenizer
from sklearn.model_selection import train_test_split

# Load dataset
df = pd.read_csv("protein_structures.csv")

# Define label mapping
label_map = {"H": 0, "B": 1, "C": 2}  # 3-class labels

# Load ProtBERT tokenizer
model_name = "Rostlab/prot_bert"
tokenizer = AutoTokenizer.from_pretrained(model_name, do_lower_case=False)

def tokenize_and_align_labels(sequence, structure):
    formatted_sequence = " ".join(sequence)
    tokenized_inputs = tokenizer(formatted_sequence, padding="max_length", truncation=True, max_length=512)

    structure_labels = [label_map.get(res, -100) for res in structure]

    # Align labels with tokenized input
    attention_mask = tokenized_inputs["attention_mask"]
    aligned_labels = [-100] * len(attention_mask)  # Initialize ignored tokens
    
    seq_index = 0
    for i, mask in enumerate(attention_mask):
        if mask == 1 and seq_index < len(structure_labels):
            aligned_labels[i] = structure_labels[seq_index]
            seq_index += 1
    
    return dict(tokenized_inputs), torch.tensor(aligned_labels)

# Process dataset
tokenized_sequences = []
label_sequences = []

for _, row in df.iterrows():
    inputs, labels = tokenize_and_align_labels(row["Primary Sequence"], row["Secondary Structure"])
    tokenized_sequences.append(inputs)  # Store as dictionary
    label_sequences.append(labels)

# Split dataset
train_inputs, test_inputs, train_labels, test_labels = train_test_split(tokenized_sequences, label_sequences, test_size=0.2, random_state=42)

# Save preprocessed data correctly (as dictionaries)
torch.save((train_inputs, train_labels), "train_data.pt")
torch.save((test_inputs, test_labels), "test_data.pt")

print("✅ Data prepared for ProtBERT fine-tuning!")
