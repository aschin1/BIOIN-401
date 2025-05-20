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

    return sov_total / len_total if len_total > 0 else 0.0

# Load model and data
model_path = "./fine_tuned_protbert"
tokenizer = AutoTokenizer.from_pretrained(model_path)
model = AutoModelForTokenClassification.from_pretrained(model_path)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

test_inputs, test_labels = torch.load("test_data.pt")

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

        for i in range(len(labels)):
            label_seq = labels[i].cpu().numpy()
            pred_seq = preds[i].cpu().numpy()
            filtered = [(p, l) for p, l in zip(pred_seq, label_seq) if l != -100]
            if filtered:
                p_clean, l_clean = zip(*filtered)
            else:
                p_clean, l_clean = [], []
            all_preds.append(p_clean)
            all_labels.append(l_clean)

# Per-sequence evaluation
sample_accuracies = []
sample_sovs = []
below_15_data = []

for i, (pred_seq, true_seq) in enumerate(zip(all_preds, all_labels)):
    length = len(true_seq)
    if length == 0:
        continue

    correct = sum(1 for p, l in zip(pred_seq, true_seq) if p == l)
    accuracy = correct / length
    sample_accuracies.append(accuracy)

    sov = compute_sov(pred_seq, true_seq)
    sample_sovs.append(sov)

    if accuracy <= 0.65:
        input_ids = test_inputs[i]["input_ids"]
        sequence = tokenizer.decode(input_ids, skip_special_tokens=True).replace(" ", "")

        pred_labels = ''.join(['H' if x == 0 else 'B' if x == 1 else 'C' for x in pred_seq])
        true_labels = ''.join(['H' if x == 0 else 'B' if x == 1 else 'C' for x in true_seq])

        min_len = min(len(sequence), len(pred_labels), len(true_labels))
        below_15_data.append([f"sample_{i}", sequence[:min_len], pred_labels[:min_len], true_labels[:min_len], min_len])

        output_dir = "low_accuracy_txt"
        os.makedirs(output_dir, exist_ok=True)
        with open(os.path.join(output_dir, f"sample_{i}.txt"), "w") as f:
            f.write(f">sample_{i}\n{sequence}\n")
            f.write(f">predicted\n{pred_labels}\n")
            f.write(f">observed\n{true_labels}\n")
        print(f"⚠️  Saved low accuracy: sample_{i}, acc={accuracy:.3f}, len={min_len}")

if below_15_data:
    with open("low_accuracy_predictions.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Accession", "Sequence", "Predicted Structure", "Actual Structure", "Length"])
        writer.writerows(below_15_data)
    print(f"💾 Saved {len(below_15_data)} low-accuracy samples to CSV")

# Summary
print(f"90-100% Accuracy: {sum(a >= 0.90 for a in sample_accuracies)}")
print(f"80-89% Accuracy: {sum(0.80 <= a < 0.90 for a in sample_accuracies)}")
print(f"70-79% Accuracy: {sum(0.70 <= a < 0.80 for a in sample_accuracies)}")
print(f"60-69% Accuracy: {sum(0.60 <= a < 0.70 for a in sample_accuracies)}")
print(f"50-59% Accuracy: {sum(0.50 <= a < 0.60 for a in sample_accuracies)}")
print(f"40-49% Accuracy: {sum(0.40 <= a < 0.50 for a in sample_accuracies)}")
print(f"30-39% Accuracy: {sum(0.30 <= a < 0.40 for a in sample_accuracies)}")
print(f"20-29% Accuracy: {sum(0.20 <= a < 0.30 for a in sample_accuracies)}")
print(f"10-19% Accuracy: {sum(0.10 <= a < 0.20 for a in sample_accuracies)}")
print(f"0-9% Accuracy: {sum(a < 0.10 for a in sample_accuracies)}")

valid_preds = [p for seq in all_preds for p in seq]
valid_labels = [l for seq in all_labels for l in seq]

print(f"\n✅ Overall Accuracy: {accuracy_score(valid_labels, valid_preds):.4f}")
print(f"🎯 Overall SOV: {compute_sov(valid_preds, valid_labels):.4f}")

print("\nClassification Report:")
print(classification_report(valid_labels, valid_preds, target_names=["H (Helix)", "B (Beta-Sheet)", "C (Coil)"], digits=4))

conf_matrix = confusion_matrix(valid_labels, valid_preds)
plt.figure(figsize=(5, 4))
sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", xticklabels=["H", "B", "C"], yticklabels=["H", "B", "C"])
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Confusion Matrix")
plt.tight_layout()
plt.show()
