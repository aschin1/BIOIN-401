def compute_structure_percent(structure_string):
    """
    Computes the percentage of H (Helix), B (Beta-Sheet), and C (Coil) in a secondary structure string.
    
    :param structure_string: A string of secondary structure labels (e.g., "HHHCCCB")
    :return: Dictionary with percentage of H, B, and C
    """
    total = len(structure_string)
    if total == 0:
        return {"H": 0.0, "B": 0.0, "C": 0.0}

    counts = {"H": 0, "B": 0, "C": 0}
    for char in structure_string:
        if char in counts:
            counts[char] += 1

    percentages = {key: (value / total) * 100 for key, value in counts.items()}
    return percentages

if __name__ == "__main__":
    structure = input("Enter the secondary structure string (e.g., 'HHHCCCB'):").strip().upper()
    percentages = compute_structure_percent(structure)

    print("\nSecondary Structure Composition:")
    print(f"Helix (H): {percentages['H']:.2f}%")
    print(f"Beta-Sheet (B): {percentages['B']:.2f}%")
    print(f"Coil (C): {percentages['C']:.2f}%")
