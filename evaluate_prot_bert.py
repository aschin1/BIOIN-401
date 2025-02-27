import os
import torch
from transformers import AutoModelForTokenClassification, AutoTokenizer
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader, Dataset

# ✅ Ensure the fine-tuned model exists
model_path = "./fine_tuned_protbert"

if not os.path.exists(model_path):
    raise FileNotFoundError(f"❌ Model directory '{model_path}' not found. Run 'fine_tune_prot_bert.py' first!")

# ✅ Load the fine-tuned model
tokenizer = AutoTokenizer.from_pretrained(model_path)
model = AutoModelForTokenClassification.from_pretrained(model_path)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

print("✅ Fine-tuned model loaded successfully!")

# ✅ Load test dataset
test_inputs, test_labels = torch.load("test_data.pt")

# Define custom dataset class
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

# Create test DataLoader
test_dataset = ProteinDataset(test_inputs, test_labels)
test_loader = DataLoader(test_dataset, batch_size=8, shuffle=False)

print(f"✅ Loaded {len(test_inputs)} test samples.")

# ✅ Evaluate the Model
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

# Remove ignored (-100) labels
valid_preds = [p for p, l in zip(all_preds, all_labels) if l != -100]
valid_labels = [l for l in all_labels if l != -100]

# ✅ Compute Accuracy
accuracy = accuracy_score(valid_labels, valid_preds)
print(f"✅ Test Accuracy: {accuracy:.4f}")

# ✅ Generate Classification Report
label_names = ["H (Helix)", "B (Beta-Sheet)", "C (Coil)"]
print("\nClassification Report:")
print(classification_report(valid_labels, valid_preds, target_names=label_names))

# ✅ Compute and Plot Confusion Matrix
conf_matrix = confusion_matrix(valid_labels, valid_preds)
plt.figure(figsize=(5, 4))
sns.heatmap(conf_matrix, annot=True, fmt="d", cmap="Blues", xticklabels=label_names, yticklabels=label_names)
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Confusion Matrix")
plt.show()
