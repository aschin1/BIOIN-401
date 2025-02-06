import subprocess
import os

def run_stride(pdb_file):
    """Runs STRIDE to analyze secondary structure from a PDB file."""
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file '{pdb_file}' not found.")
        return None

    command = ["stride", pdb_file]
    try:
        output = subprocess.check_output(command, text=True)
        if not output.strip():
            print("Error: STRIDE did not return any output. Check if the PDB file is valid.")
            return None
        return output
    except subprocess.CalledProcessError as e:
        print(f"Error running STRIDE: {e}")
        return None

def parse_stride_output(output):
    """Parses STRIDE output to extract residue sequence and secondary structure labels."""
    residues = []  # Stores amino acid residues
    sec_structure = []  # Stores corresponding secondary structure labels

    lines = output.split("\n")

    for line in lines:
        if line.startswith("ASG"):  # STRIDE uses "ASG" for secondary structure assignment
            parts = line.split()
            if len(parts) >= 6:
                try:
                    residue = parts[1][0]  # Column 2 = Amino Acid Code
                    structure_label = parts[5]  # Column 6 = Secondary structure type
                    
                    residues.append(residue)  # Add amino acid to sequence
                    sec_structure.append(structure_label)  # Add structure type
                except ValueError:
                    continue

    if not residues or not sec_structure:
        print("Error: No secondary structure information extracted from STRIDE output.")
        return None, None

    return residues, sec_structure

def display_horizontal_output(residues, sec_structure, block_size=50):
    """Formats and prints the sequence and secondary structure with perfect alignment."""
    print("\nPrimary Sequence and Secondary Structure (Horizontal Format):\n")

    for i in range(0, len(residues), block_size):
        chunk_residues = residues[i:i+block_size]
        chunk_structure = sec_structure[i:i+block_size]
        
        # Create sequence and structure lines with exact alignment
        seq_line = f"SEQ  {i+1:<5} " + "".join(chunk_residues) + f"  {min(i+block_size, len(residues))}"
        str_line = f"STR  {i+1:<5} " + "".join(chunk_structure)  # Match spaces to seq_line

        print(seq_line)
        print(str_line)
        print()  # Blank line between blocks

if __name__ == "__main__":
    pdb_path = input("Enter the full path to the PDB file: ").strip()

    stride_output = run_stride(pdb_path)
    
    if stride_output:
        residues, sec_structure = parse_stride_output(stride_output)

        if residues and sec_structure:
            display_horizontal_output(residues, sec_structure)
        else:
            print("Error: No secondary structure data extracted.")
    else:
        print("Error: STRIDE execution failed.")
