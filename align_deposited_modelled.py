from Bio import pairwise2

# Step 1: Enter sequences manually
vadar_seq = input("Enter VADAR modeled sequence (e.g., MSIVSY...): ").strip().upper()
deposited_seq = input("Enter deposited sequence (same as FASTA): ").strip().upper()
model_preds = input("Enter model predictions (1 char per residue, same length as deposited): ").strip().upper()

if len(deposited_seq) != len(model_preds):
    print("❌ Error: Deposited sequence and model prediction lengths must match.")
    exit()

# Step 2: Align sequences
alignment = pairwise2.align.globalxx(vadar_seq, deposited_seq, one_alignment_only=True)[0]

vadar_aligned = alignment.seqA
deposited_aligned = alignment.seqB

# Step 3: Print aligned sequences
print("\nAligned Sequences:")
print("VADAR     :", vadar_aligned)
print("Deposited :", deposited_aligned)

# Step 4: Find match positions in deposited sequence
deposited_index = 0
matched_indices = []

for v_char, d_char in zip(vadar_aligned, deposited_aligned):
    if d_char != "-":
        if v_char == d_char:
            matched_indices.append(deposited_index)
        deposited_index += 1

# Step 5: Extract matching model predictions
model_aligned = "".join([model_preds[i] for i in matched_indices])

# Step 6: Show output
print("\nMatched Residue Indices in Deposited Sequence:", matched_indices)
print("Extracted Model Predictions at Those Positions:")
print(model_aligned)
