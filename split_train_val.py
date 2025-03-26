import torch
from sklearn.model_selection import train_test_split

# Load full dataset
full_inputs, full_labels = torch.load("train_data.pt")

# Split data (80% training, 20% validation)
train_inputs, val_inputs, train_labels, val_labels = train_test_split(
    full_inputs, full_labels, test_size=0.2, random_state=42
)

# Save training and validation datasets separately
torch.save((train_inputs, train_labels), "train.pt")
torch.save((val_inputs, val_labels), "validation.pt")

print(f"✅ Training set: {len(train_inputs)} samples")
print(f"✅ Validation set: {len(val_inputs)} samples")
print("✅ Data successfully split and saved!")
