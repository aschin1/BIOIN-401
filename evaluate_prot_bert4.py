# save_evaluation_results.py

import os
import torch
import numpy as np
import csv
from transformers import AutoModelForTokenClassification, AutoTokenizer
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader, Dataset
from collections import defaultdict

def compute_sov(preds, labels, label_set={0, 1, 2}):
    def get_segments(sequence):
        segments = []
        if not sequence:
            return segments
        start = 0
        current_label = sequence[0]
        for i in range(1, len(sequence)):
            if sequence[i] != current_label:
                segments.append((start, i - 1, current_label))
                start = i
                current_label = sequence[i]
        segments.append((start, len(sequence) - 1, current_label))
        return segments

    true_segments = get_segments(labels)
    pred_segments = get_segments(preds)

    sov_total = 0
    len_total = 0

    for label in label_set:
        true_segs = [seg for seg in true_segments if seg[2] == label]
        pred_segs = [seg for seg in pred_segments if seg[2] == label]

        for t_start, t_end, _ in true_segs:
            t_len = t_end - t_start + 1
            len_total += t_len

            overlaps = []
            for p_start, p_end, _ in pred_segs:
                ov_start = max(t_start, p_start)
                ov_end = min(t_end, p_end)
                if ov_start <= ov_end:
                    ov_len = ov_end - ov_start + 1
                    union_len = max(t_end, p_end) - min(t_start, p_start) + 1
                    delta = min([
                        ov_len,
                        t_end - t_start + 1 - ov_len,
                        p_end - p_start + 1 - ov_len
                    ])
                    sov_contrib = ((ov_len + delta) / union_len) * t_len
                    overlaps.append(sov_contrib)

            if overlaps:
                sov_total += max(overlaps)
            else:
                sov_total += 0

    return sov_total / len_total if len_total > 0 else 0.0

# Load fine-tuned model
model_path = "./fine_tuned_protbert"
tokenizer = AutoTokenizer.from_pretrained(model_path)
model = AutoModelForTokenClassification.from_pretrained(model_path)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

# Load test data
test_inputs, test_labels = torch.load("test_data.pt")

# Dataset wrapper
class ProteinDataset(Dataset):
    def __init__(self, inputs, labels):
        self.inputs = inputs
        self.labels = labels

    def __len__(self):
        return len(self.inputs)

    def __getitem__(self, idx):
        inputs = {key: torch.tensor(val).squeeze() for key, val in self.inputs[idx].items()}
        labels = self.labels[idx]
        return inputs, labels

test_dataset = ProteinDataset(test_inputs, test_labels)
test_loader = DataLoader(test_dataset, batch_size=8, shuffle=False)

# Predict
all_preds = []
all_labels = []

model.eval()
with torch.no_grad():
    for batch in test_loader:
        inputs, labels = batch
        inputs = {key: val.to(device) for key, val in inputs.items()}
        labels = labels.to(device)
        outputs = model(**inputs).logits
        _, preds = torch.max(outputs, 2)
        all_preds.extend(preds.cpu().numpy().flatten())
        all_labels.extend(labels.cpu().numpy().flatten())

# Filter padding
valid_preds = [p for p, l in zip(all_preds, all_labels) if l != -100]
valid_labels = [l for l in all_labels if l != -100]

# Evaluate and collect poor-performing examples
sample_accuracies = []
sample_sovs = []
below_60_data = []

start_idx = 0
for i, label_seq in enumerate(test_labels):
    length = len(label_seq)
    end_idx = start_idx + length

    true_seq = valid_labels[start_idx:end_idx]
    pred_seq = valid_preds[start_idx:end_idx]

    correct = sum(1 for p, l in zip(pred_seq, true_seq) if p == l)
    accuracy = correct / length if length > 0 else 0
    sample_accuracies.append(accuracy)

    sov = compute_sov(pred_seq, true_seq)
    sample_sovs.append(sov)

    if accuracy <= 0.60:
        input_ids = test_inputs[i]["input_ids"]
        tokens = tokenizer.convert_ids_to_tokens(input_ids)
        amino_acids = [t.replace("▁", "") for t in tokens if t not in ("[CLS]", "[SEP]", "-", "<pad>")]
        sequence = ''.join(amino_acids)

        pred_labels = ''.join(['H' if x == 0 else 'B' if x == 1 else 'C' for x in pred_seq])
        true_labels = ''.join(['H' if x == 0 else 'B' if x == 1 else 'C' for x in true_seq])
        accession = f"sample_{i}"

        below_60_data.append([accession, sequence, pred_labels, true_labels])

    start_idx = end_idx

# Save poor predictions
if below_60_data:
    with open("low_accuracy_predictions.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Accession", "Sequence", "Predicted Structure", "Actual Structure"])
        writer.writerows(below_60_data)
    print(f"💾 Saved {len(below_60_data)} samples <60% accuracy to 'low_accuracy_predictions.csv'")

# Summary
print(f"\n📊 >=90% Accuracy: {sum(a >= 0.90 for a in sample_accuracies)}")
print(f"✅ Fully Accurate: {sum(a == 1.0 for a in sample_accuracies)}")
print(f" 80-89% Accuracy: {sum(80 <= a*100 < 90 for a in sample_accuracies)}")
print(f" 70-79% Accuracy: {sum(70 <= a*100 < 80 for a in sample_accuracies)}")
print(f" 60-69% Accuracy: {sum(60 <= a*100 < 70 for a in sample_accuracies)}")
print(f" 50-59% Accuracy: {sum(50 <= a*100 < 60 for a in sample_accuracies)}")
print(f" 40-49% Accuracy: {sum(40 <= a*100 < 50 for a in sample_accuracies)}")
print(f" 30-39% Accuracy: {sum(30 <= a*100 < 40 for a in sample_accuracies)}")
print(f" 20-29% Accuracy: {sum(20 <= a*100 < 30 for a in sample_accuracies)}")
print(f" 10-19% Accuracy: {sum(10 <= a*100 < 20 for a in sample_accuracies)}")
print(f" 0-9% Accuracy: {sum(a*100 < 10 for a in sample_accuracies)}")



accuracy = accuracy_score(valid_labels, valid_preds)
sov_score = compute_sov(valid_preds, valid_labels)
print(f"\n✅ Overall Test Accuracy: {accuracy:.4f}")
print(f"🎯 Overall SOV Score: {sov_score:.4f}")

# Classification report and confusion matrix
label_names = ["H (Helix)", "B (Beta-Sheet)", "C (Coil)"]
print("\nClassification Report:")
print(classification_report(valid_labels, valid_preds, target_names=label_names, digits=4))

conf_matrix = confusion_matrix(valid_labels, valid_preds)
plt.figure(figsize=(5, 4))
sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", xticklabels=label_names, yticklabels=label_names)
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Confusion Matrix")
plt.tight_layout()
plt.show()
