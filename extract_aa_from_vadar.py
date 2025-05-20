import subprocess
import sys
import os
import re

def extract_sequence_from_vadar_output(vadar_text):
    lines = vadar_text.strip().splitlines()
    sequence = ""
    structure = ""

    for line in lines:
        if line.startswith("SEQ"):
            match = re.match(r"SEQ\s+\d+\s+([A-Z]+)\s+\d+", line)
            if match:
                sequence += match.group(1)
            else:
                parts = line.split()
                if len(parts) >= 2:
                    sequence += parts[2] if len(parts) > 2 else parts[1]
        elif line.startswith("STR"):
            match = re.match(r"STR\s+\d+\s+([A-Z]+)\s+\d+", line)
            if match:
                structure += match.group(1)
            else:
                parts = line.split()
                if len(parts) >= 2:
                    structure += parts[2] if len(parts) > 2 else parts[1]

    return sequence, structure


def run_get_secondary_struct(pdb_path):
    try:
        result = subprocess.run(
            ["python", "get_secondary_struct_dssp.py"],
            input=pdb_path + "\n",
            text=True,
            capture_output=True
        )
        if result.returncode != 0:
            print("❌ Error running get_secondary_struct_dssp.py:", result.stderr)
            return None
        return result.stdout
    except Exception as e:
        print("❌ Exception occurred:", str(e))
        return None

def main():
    if len(sys.argv) != 2:
        print("Usage: python extract_vadar_sequence.py <pdb_file>")
        sys.exit(1)

    pdb_path = sys.argv[1]
    if not os.path.exists(pdb_path):
        print(f"❌ File not found: {pdb_path}")
        sys.exit(1)

    vadar_output = run_get_secondary_struct(pdb_path)
    if not vadar_output:
        sys.exit(1)

    sequence, structure = extract_sequence_from_vadar_output(vadar_output)
    print("\n✅ Extracted amino acid sequence from VADAR:")
    print(sequence)
    print("\n✅ Extracted secondary structure from VADAR:")
    print(structure)

if __name__ == "__main__":
    main()

