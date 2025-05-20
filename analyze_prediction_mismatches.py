def analyze_prediction_file(txt_path):
    with open(txt_path, 'r') as file:
        lines = [line.strip() for line in file.readlines() if line.strip() and not line.startswith('>')]

    if len(lines) != 3:
        print("❌ File format error: expected 3 lines (sequence, predicted, observed)")
        return

    sequence, predicted, observed = lines

    if not (len(predicted) == len(observed) == len(sequence)):
        print("❌ Length mismatch between sequence and predictions")
        return

    mismatches = []
    for i, (pred, obs) in enumerate(zip(predicted, observed)):
        if pred != obs:
            mismatches.append((i + 1, sequence[i], pred, obs))  # 1-based index

    print("\n🧪 Mismatch Report:")
    if not mismatches:
        print("✅ No mismatches detected!")
    else:
        print(f"Found {len(mismatches)} mismatches:")
        print(f"{'Pos':>4}  {'AA':>2}  {'Pred':>5}  {'Obs':>5}")
        print("-" * 25)
        for idx, aa, pred, obs in mismatches:
            print(f"{idx:>4}  {aa:>2}  {pred:>5}  {obs:>5}")

        # Group into contiguous mismatch regions
        print("\n📦 Contiguous Mismatch Regions:")
        regions = []
        current_region = []
        for i, mismatch in enumerate(mismatches):
            if not current_region:
                current_region.append(mismatch)
            else:
                if mismatch[0] == current_region[-1][0] + 1:
                    current_region.append(mismatch)
                else:
                    regions.append(current_region)
                    current_region = [mismatch]
        if current_region:
            regions.append(current_region)

        for region in regions:
            start = region[0][0]
            end = region[-1][0]
            length = end - start + 1
            print(f"Positions {start}-{end} (length {length})")

        # Write marked alignment to file
        def write_alignment_with_mismatches(sequence, predicted, observed, output_file):
            marker_line = ['^' if p != o else ' ' for p, o in zip(predicted, observed)]
            with open(output_file, 'w') as f:
                f.write("SEQ : " + sequence + "\n")
                f.write("PRED: " + predicted + "\n")
                f.write("OBS : " + observed + "\n")
                f.write("      " + ''.join(marker_line) + "\n")
            print(f"\n📄 Alignment with mismatches saved to: {output_file}")

        output_file = txt_path.replace(".txt", "_alignment.txt")
        write_alignment_with_mismatches(sequence, predicted, observed, output_file)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python analyze_prediction_mismatches.py <file.txt>")
    else:
        analyze_prediction_file(sys.argv[1])
