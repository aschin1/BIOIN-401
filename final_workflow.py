import os
import torch
import pandas as pd
from transformers import AutoModelForTokenClassification, AutoTokenizer
from Bio import SeqIO

# ✅ Load the fine-tuned ProtBERT model
MODEL_PATH = "./fine_tuned_protbert"

if not os.path.exists(MODEL_PATH):
    raise FileNotFoundError(f"❌ Model directory '{MODEL_PATH}' not found. Train the model first!")

tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH)
model = AutoModelForTokenClassification.from_pretrained(MODEL_PATH)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

print("✅ Fine-tuned model loaded successfully!")

# ✅ Load existing dataset (if available)
DATASET_PATH = "protein_structures_50.csv"
if os.path.exists(DATASET_PATH):
    df = pd.read_csv(DATASET_PATH)
else:
    df = None

# ✅ Function to retrieve secondary structure from dataset
def get_secondary_structure_from_dataset(sequence):
    if df is not None:
        match = df[df["Primary Sequence"] == sequence]
        if not match.empty:
            return match["Secondary Structure"].values[0]
    return None

# ✅ Function to tokenize input sequence and predict secondary structure
def predict_secondary_structure(sequence, model, tokenizer, max_len=512):
    """
    Predicts the secondary structure of a given protein sequence using ProtBERT.

    :param sequence: Protein sequence (string).
    :param model: Fine-tuned ProtBERT model.
    :param tokenizer: Corresponding tokenizer.
    :param max_len: Maximum token length.
    :return: Predicted secondary structure as a string.
    """
    model.eval()

    # ✅ Fix: Ensure correct tokenization without padding errors
    tokenized_inputs = tokenizer(" ".join(sequence),
                                 padding="max_length",
                                 truncation=True,
                                 max_length=max_len,
                                 return_tensors="pt")

    input_ids = tokenized_inputs["input_ids"].to(device)
    attention_mask = tokenized_inputs["attention_mask"].to(device)

    with torch.no_grad():
        logits = model(input_ids, attention_mask=attention_mask).logits
        predictions = torch.argmax(logits, dim=2).cpu().numpy().flatten()

    # ✅ Fix: Trim predictions to match sequence length (ignore special tokens)
    valid_predictions = predictions[:len(sequence)]  # Ensures correct length

    # Convert numeric labels back to H/B/C characters
    label_map = {0: "H", 1: "B", 2: "C"}
    predicted_structure = "".join([label_map[p] for p in valid_predictions])

    return predicted_structure


# ✅ Function to extract sequence from a FASTA file
def extract_sequence_from_fasta(fasta_file):
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            return str(record.seq)
    return None

# ✅ Function to display aligned sequences
def display_aligned_sequences(primary_seq, secondary_struct):
    print("\nPrimary Sequence and Predicted Secondary Structure:")
    print("SEQ  :", primary_seq)
    print("STR  :", secondary_struct)
    print("\n" + "-" * 60)

# ✅ User Interaction
def main():
    print("\nProtein Secondary Structure Prediction")
    print("1. Enter a protein sequence")
    print("2. Upload a FASTA file")
    choice = input("\nChoose an option (1 or 2): ").strip()

    if choice == "1":
        sequence = input("\nEnter your protein sequence: ").strip()
    elif choice == "2":
        fasta_file = input("\nEnter the path to your FASTA file: ").strip()
        if not os.path.exists(fasta_file):
            print("❌ Error: File not found!")
            return
        sequence = extract_sequence_from_fasta(fasta_file)
        print(f"\n✅ Extracted sequence: {sequence}")
    else:
        print("❌ Invalid choice. Exiting.")
        return

    # ✅ Initialize `existing_structure` to prevent errors
    existing_structure = get_secondary_structure_from_dataset(sequence)

    if existing_structure:
        print("\n✅ Secondary structure found in dataset!")
    else:
        print("\n🔍 Secondary structure not found. Predicting using ProtBERT...")
        existing_structure = predict_secondary_structure(sequence, model, tokenizer)

    # ✅ Display aligned sequences
    display_aligned_sequences(sequence, existing_structure)

if __name__ == "__main__":
    main()
