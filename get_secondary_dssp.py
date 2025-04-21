import os
import subprocess

# DSSP 8-class to 3-class mapping
DSSP_TO_3CLASS = {
    'H': 'H',  # Alpha-helix → Helix (H)
    'G': 'H',  # 3_10-helix → Helix (H)
    'I': 'H',  # Pi-helix → Helix (H)
    'E': 'B',  # Beta-sheet → Beta-strand (B)
    'B': 'B',  # Beta-bridge → Beta-strand (B)
    'T': 'C',  # Turn → Coil (C)
    'S': 'C',  # Bend → Coil (C)
    '-': 'C'   # Loop/Irregular → Coil (C)
}

def run_dssp(pdb_file, dssp_executable="mkdssp", output_file="dssp_output.txt"):
    """Runs DSSP as an external command and saves output to a file."""
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file '{pdb_file}' not found.")
        return None
    
    try:
        with open(output_file, "w") as out:
            subprocess.run([dssp_executable, "-i", pdb_file, "-o", output_file], stdout=out, text=True, check=True)
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"Error running DSSP: {e}")
        return None

def parse_dssp_output(output_file):
    """Parses DSSP output file and extracts residues and secondary structure labels."""
    if not os.path.exists(output_file):
        print(f"Error: DSSP output file '{output_file}' not found.")
        return None, None
    
    residues = []
    sec_structure = []
    
    try:
        with open(output_file, "r") as f:
            lines = f.readlines()
        
        data_start = False
        for line in lines:
            if line.startswith("  #  RESIDUE"):  # DSSP header line before data
                data_start = True
                continue
            
            if not data_start or len(line) < 17:
                continue
            
            residue = line[13]  # Amino acid one-letter code
            structure_label = line[16]  # Secondary structure DSSP code
            
            # Convert DSSP 8-class to 3-class
            structure_label = DSSP_TO_3CLASS.get(structure_label, 'C')
            
            residues.append(residue)
            sec_structure.append(structure_label)
        
        if not residues or not sec_structure:
            print("Error: No secondary structure data extracted from DSSP output.")
            return None, None
        
        return residues, sec_structure
    except Exception as e:
        print(f"Error parsing DSSP output: {e}")
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
    
    dssp_output = run_dssp(pdb_path)
    if not dssp_output:
        print("Error: DSSP execution failed.")
        exit(1)
    
    residues, sec_structure = parse_dssp_output(dssp_output)
    if not residues or not sec_structure:
        print("Error: No secondary structure data extracted.")
        exit(1)
    
    display_horizontal_output(residues, sec_structure)
