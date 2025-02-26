import os
import subprocess

# Mapping of three-letter amino acid codes to one-letter codes
amino_acid_map = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", 
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", 
    "TYR": "Y", "VAL": "V"
}

def run_phipsi(pdb_file, phipsi_executable="./phipsi", output_file="phipsi_output.txt"):
    """Runs phipsi as an external command and saves output to a file."""
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file '{pdb_file}' not found.")
        return None

    if not os.path.exists(phipsi_executable):
        print(f"Error: '{phipsi_executable}' executable not found. Ensure it has execution permissions.")
        return None

    try:
        with open(output_file, "w") as out:
            subprocess.run([phipsi_executable, pdb_file], stdout=out, text=True, check=True)
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"Error running phipsi: {e}")
        return None

def run_define2(phipsi_output_file, define2_executable="./define2", output_file="secondary_structure.txt"):
    """Runs define2 on the phipsi output file to classify secondary structures."""
    if not os.path.exists(phipsi_output_file):
        print(f"Error: Phipsi output file '{phipsi_output_file}' not found.")
        return None

    if not os.path.exists(define2_executable):
        print(f"Error: '{define2_executable}' executable not found. Ensure it has execution permissions.")
        return None

    try:
        with open(output_file, "w") as out:
            subprocess.run([define2_executable, "-f", phipsi_output_file], stdout=out, text=True, check=True)
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"Error running define2: {e}")
        return None

def parse_define2_output(output_file):
    """Parses define2 output file and extracts residues and secondary structure labels."""
    if not os.path.exists(output_file):
        print(f"Error: Secondary structure file '{output_file}' not found.")
        return None, None

    residues = []
    sec_structure = []
    
    try:
        with open(output_file, "r") as f:
            lines = f.readlines()

        print("\nDEBUG: Raw define2 Output:\n", "".join(lines))  # Print file content for debugging

        separator_count = 0
        data_start = False
        for line in lines:
            line = line.strip()

            # Detect the second separator line before starting data extraction
            if line.startswith("------------------------------------"):
                separator_count += 1
                if separator_count == 2:
                    data_start = True
                continue

            if not data_start or line == "":
                continue  # Ignore lines before the structured data starts

            # Process structured data (Residue Number, Name, Secondary Structure)
            parts = line.split()
            
            # Ensure the line has at least three valid columns (NUM, NAME, STRUCTURE)
            if len(parts) >= 3 and parts[0].isdigit():
                residue = amino_acid_map.get(parts[1].strip().upper(), "X")  # Convert 3-letter to 1-letter, default to X if not found
                structure_label = parts[2].strip()  # Secondary Structure (Column 3)

                residues.append(residue)
                sec_structure.append(structure_label)

        if not residues or not sec_structure:
            print("Error: No valid secondary structure data extracted from define2 output.")
            return None, None

        return residues, sec_structure
    except Exception as e:
        print(f"Error parsing define2 output: {e}")
        return None, None

def display_horizontal_output(residues, sec_structure, block_size=50):
    """Formats and prints the sequence and secondary structure horizontally."""
    print("\nPrimary Sequence and Secondary Structure (Horizontal Format):\n")

    seq_label = "SEQ  "
    str_label = "STR  "

    for i in range(0, len(residues), block_size):
        chunk_residues = residues[i:i+block_size]
        chunk_structure = sec_structure[i:i+block_size]

        formatted_residues = "".join(chunk_residues)
        seq_line = f"{seq_label}{i+1:<5} " + formatted_residues + f"  {min(i+block_size, len(residues))}"
        
        formatted_structure = "".join(chunk_structure)
        str_line = f"{str_label}{' ' * 6}" + formatted_structure

        print(seq_line)
        print(str_line)
        print()  # Blank line between blocks

if __name__ == "__main__":
    pdb_path = input("Enter the full path to the PDB file: ").strip()

    phipsi_output = run_phipsi(pdb_path)
    if not phipsi_output:
        print("Error: Phipsi execution failed.")
        exit(1)
    
    define2_output = run_define2(phipsi_output)
    if not define2_output:
        print("Error: Define2 execution failed.")
        exit(1)
    
    residues, sec_structure = parse_define2_output(define2_output)
    if not residues or not sec_structure:
        print("Error: No secondary structure data extracted.")
        exit(1)
    
    display_horizontal_output(residues, sec_structure)
