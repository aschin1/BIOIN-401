import pandas as pd

# Load the dataset
csv_path = "/mnt/data/protein_structures_rand_50.csv"
df = pd.read_csv(csv_path)

# Ensure dataset has the required column
if "Primary Sequence" not in df.columns:
    raise ValueError("The dataset does not contain a 'Primary Sequence' column.")

# Define train-test split ratio (assumed 80-20 split)
train_size = int(0.8 * len(df))

# Extract training and test sequences
train_sequences = set(df.iloc[:train_size]["Primary Sequence"])
test_sequences = set(df.iloc[train_size:]["Primary Sequence"])

# Find overlapping sequences
overlapping_sequences = train_sequences.intersection(test_sequences)

# Calculate overlap percentage
overlap_count = len(overlapping_sequences)
total_test_samples = len(test_sequences)
overlap_percentage = (overlap_count / total_test_samples) * 100 if total_test_samples > 0 else 0

# Display results
import ace_tools as tools

df_results = pd.DataFrame({
    "Metric": ["Total Training Sequences", "Total Test Sequences", "Overlapping Sequences", "Overlap Percentage"],
    "Value": [len(train_sequences), len(test_sequences), overlap_count, f"{overlap_percentage:.2f}%"]
})

tools.display_dataframe_to_user(name="Train-Test Overlap Check", dataframe=df_results)
