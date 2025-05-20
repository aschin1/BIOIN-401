def compute_similarity(seq1, seq2):
    """
    Computes the percentage of identical characters between two sequences.
    Sequences must be of equal length.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length to compute similarity.")

    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    similarity = matches / len(seq1) * 100
    return similarity

if __name__ == "__main__":
    seq1 = input("Enter the first sequence: ").strip()
    seq2 = input("Enter the second sequence: ").strip()

    try:
        similarity = compute_similarity(seq1, seq2)
        print(f"Sequence Similarity: {similarity:.2f}%")
    except ValueError as e:
        print(f"Error: {e}")
