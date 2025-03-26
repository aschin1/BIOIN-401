import os
import torch
import pandas as pd
from transformers import AutoModelForTokenClassification, AutoTokenizer
from Bio import SeqIO
from accession_retrieval4 import get_fasta_metadata
from get_alphafold_data4 import fetch_alphafold_pdb, save_pdb_file
from get_secondary_struct4 import run_phipsi, run_define2, parse_define2_output

# ✅ Load the fine-tuned ProtBERT model
MODEL_PATH = "./fine_tuned_protbert"

if not os.path.exists(MODEL_PATH):
    raise FileNotFoundError(f"❌ Model directory '{MODEL_PATH}' not found. Train the model first!")

tokenizer = AutoTokenizer.from_pretrained(MODEL_PATH)
model = AutoModelForTokenClassification.from_pretrained(MODEL_PATH)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

print("✅ Fine-tuned model loaded successfully!")

def sliding_window(sequence, window_size=512, num_segments=10):
    """
    Generate overlapping segments of a protein sequence.
    """
    seq_length = len(sequence)
    if seq_length < window_size:
        return [sequence]  # If sequence is too short, return as a single segment

    # Compute step size for even overlap
    step_size = max(1, (seq_length - window_size) // (num_segments - 1))
    segments = [sequence[i : i + window_size] for i in range(0, seq_length - window_size + 1, step_size)]
    return segments[:num_segments]  # Ensure exactly `num_segments` windows

def sliding_window_prediction(sequence, model, tokenizer, window_size=512, num_segments=10):
    """
    Apply sliding window to predict secondary structure for long sequences.
    """
    label_map = {0: "H", 1: "B", 2: "C"}
    segments = sliding_window(sequence, window_size, num_segments)
    predictions = []

    for segment in segments:
        tokenized_inputs = tokenizer(" ".join(segment),
                                     padding="max_length",
                                     truncation=True,
                                     max_length=window_size,
                                     return_tensors="pt")
        input_ids = tokenized_inputs["input_ids"].to(device)
        attention_mask = tokenized_inputs["attention_mask"].to(device)

        with torch.no_grad():
            logits = model(input_ids, attention_mask=attention_mask).logits
            pred_labels = torch.argmax(logits, dim=2).cpu().numpy().flatten()
        
        predictions.append("".join([label_map[p] for p in pred_labels[:len(segment)]]))

    if not predictions:
        return "C" * len(sequence)  # Default to coil structure if no predictions are made

    # Merge overlapping predictions
    merged_prediction = predictions[0]
    for i in range(1, len(predictions)):
        merged_prediction += predictions[i][-(window_size // 2):]  # Overlapping merge

    return merged_prediction[:len(sequence)]

def extract_sequence_from_fasta(fasta_file):
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            return str(record.seq)
    return None

def fetch_alphafold_secondary_structure(accession):
    pdb_data = fetch_alphafold_pdb(accession)
    if not pdb_data:
        return None
    
    pdb_filename = f"pdb_files/{accession}.pdb"
    save_pdb_file(accession, pdb_data)
    
    phipsi_output = run_phipsi(pdb_filename)
    if not phipsi_output:
        return None
    
    define2_output = run_define2(phipsi_output)
    if not define2_output:
        return None
    
    residues, sec_structure = parse_define2_output(define2_output)
    return "".join(sec_structure) if sec_structure else None

def display_aligned_sequences(primary_seq, secondary_struct):
    print("\nPrimary Sequence and Predicted Secondary Structure:")
    print("SEQ  :", primary_seq)
    print("STR  :", secondary_struct)
    print("\n" + "-" * 60)

def main():
    print("\nProtein Secondary Structure Prediction")
    print("1. Enter a protein sequence")
    print("2. Upload a FASTA file")
    choice = input("\nChoose an option (1 or 2): ").strip()

    if choice == "1":
        sequence = input("\nEnter your protein sequence: ").strip()
        accession = None
    elif choice == "2":
        fasta_file = input("\nEnter the path to your FASTA file: ").strip()
        if not os.path.exists(fasta_file):
            print("❌ Error: File not found!")
            return
        accession, _ = get_fasta_metadata(fasta_file)
        sequence = extract_sequence_from_fasta(fasta_file)
        print(f"\n✅ Extracted accession number: {accession}")
    else:
        print("❌ Invalid choice. Exiting.")
        return

    print("\n🔍 Checking AlphaFold for secondary structure...")
    if accession:
        alphafold_structure = fetch_alphafold_secondary_structure(accession)
    else:
        alphafold_structure = None

    if alphafold_structure:
        print("\n✅ Retrieved secondary structure from AlphaFold!")
        predicted_structure = alphafold_structure
    else:
        print("\n🔍 AlphaFold data unavailable. Predicting using sliding window approach...")
        predicted_structure = sliding_window_prediction(sequence, model, tokenizer)

    display_aligned_sequences(sequence, predicted_structure)

if __name__ == "__main__":
    main()

