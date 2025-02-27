import torch
import torch.nn as nn
import torch.optim as optim
from transformers import AutoModelForTokenClassification, AutoTokenizer
from torch.utils.data import DataLoader, Dataset

# ✅ Load preprocessed dataset
train_inputs, train_labels = torch.load("train_data.pt")

# Define a dataset class
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

# Create DataLoader
train_dataset = ProteinDataset(train_inputs, train_labels)
train_loader = DataLoader(train_dataset, batch_size=2, shuffle=True)

print(f"✅ Loaded {len(train_inputs)} training samples.")

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

print("✅ Fine-tuning complete!")

# ✅ Save the fine-tuned model
model.save_pretrained("fine_tuned_protbert")
tokenizer.save_pretrained("fine_tuned_protbert")

print("✅ Fine-tuned model saved successfully!")
