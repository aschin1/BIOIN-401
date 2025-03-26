import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from transformers import AutoModelForTokenClassification, AutoTokenizer
from sklearn.metrics import accuracy_score
from torch.utils.data import DataLoader, Dataset
from torch.nn.utils.rnn import pad_sequence

# ✅ Load preprocessed dataset (with sliding window applied)
train_inputs, train_labels = torch.load("train.pt")
val_inputs, val_labels = torch.load("validation.pt")  # Load validation data

# Define dataset class
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

# Create DataLoader with correct batching
train_dataset = ProteinDataset(train_inputs, train_labels)
train_loader = DataLoader(train_dataset, batch_size=8, shuffle=True)
num_residues = sum(len(label) for label in train_labels)
val_dataset = ProteinDataset(val_inputs, val_labels)
val_loader = DataLoader(val_dataset, batch_size=8, shuffle=False)


print(f"✅ Loaded {len(train_inputs)} training samples and {num_residues} labels (with sliding windows).")

# ✅ Load ProtBERT model for fine-tuning
model_name = "Rostlab/prot_bert"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForTokenClassification.from_pretrained(model_name, num_labels=3)  # 3 classes: H, B, C
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

# Define loss function and optimizer
criterion = nn.CrossEntropyLoss(ignore_index=-100)  # Ignores padding tokens
optimizer = optim.AdamW(model.parameters(), lr=2e-5)

# ✅ Fine-tune the model
num_epochs = 3  # Increase for better results

for epoch in range(num_epochs):
    model.train()
    total_loss = 0

    for batch in train_loader:
        inputs, labels = batch
        inputs = {key: val.to(device) for key, val in inputs.items()}
        labels = labels.to(device)

        optimizer.zero_grad()
        outputs = model(**inputs).logits
        loss = criterion(outputs.view(-1, 3), labels.view(-1))

        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {total_loss / len(train_loader):.4f}")
    
    
    model.eval()
    val_preds = []
    val_labels_flat = []
    
    with torch.no_grad():
        for batch in val_loader:
            inputs, labels = batch
            inputs = {key: val.to(device) for key, val in inputs.items()}
            labels = labels.to(device)

            outputs = model(**inputs).logits
            _, preds = torch.max(outputs, 2)
            
            val_preds.extend(preds.cpu().numpy().flatten())
            val_labels_flat.extend(labels.cpu().numpy().flatten())

    # Remove ignored (-100) labels
    # Remove ignored (-100) labels
    valid_preds = [p for p, l in zip(val_preds, val_labels_flat) if l != -100]
    valid_labels = [l for l in val_labels_flat if l != -100]

    # Debugging output
    print(f"Validation Samples: {len(val_labels_flat)}, Filtered Samples: {len(valid_labels)}")

    if len(valid_labels) == 0 or len(valid_preds) == 0:
        print("⚠️ Warning: No valid validation labels found. Skipping accuracy computation.")
    else:
        val_accuracy = accuracy_score(valid_labels, valid_preds)

    print(f"Epoch [{epoch+1}/{num_epochs}], Validation Accuracy: {val_accuracy:.4f}")

print("✅ Fine-tuning complete!")

# ✅ Save the fine-tuned model
model.save_pretrained("fine_tuned_protbert")
tokenizer.save_pretrained("fine_tuned_protbert")

print("✅ Fine-tuned model saved successfully!")
