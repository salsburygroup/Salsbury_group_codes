import sys
import os
import subprocess
import tempfile

def read_clean_protein(protein_pdb):
    """Reads a PDB file and removes the final END statement, if present."""
    with open(protein_pdb, 'r') as f:
        lines = f.readlines()
    return [line for line in lines if not line.strip().startswith("END")]

def convert_ligand_model_to_pdb(ligand_lines, temp_pdbqt_path):
    """Write ligand PDBQT model to a temp file, fix it using Open Babel, return fixed lines."""
    with open(temp_pdbqt_path, "w") as temp_file:
        temp_file.writelines(ligand_lines)

    fixed_pdb_path = temp_pdbqt_path.replace(".pdbqt", "_fixed.pdb")

    # Run obabel to convert/fix the file
    subprocess.run([
        "obabel", "-ipdbqt", temp_pdbqt_path, "-opdb", "-O", fixed_pdb_path
    ], check=True)

    with open(fixed_pdb_path, "r") as f:
        fixed_lines = f.readlines()

    # Cleanup temp files
    os.remove(temp_pdbqt_path)
    os.remove(fixed_pdb_path)

    return fixed_lines

def split_and_combine_models(pdbqt_file, protein_pdb, output_base="complex"):
    with open(pdbqt_file, 'r') as f:
        lines = f.readlines()

    model = []
    model_index = 0
    writing_model = False

    protein_lines = read_clean_protein(protein_pdb)

    for line in lines:
        if line.startswith("MODEL"):
            writing_model = True
            model = [line]
        elif line.startswith("ENDMDL"):
            model.append(line)
            model_index += 1
            # Remove docking-specific lines (ROOT, BRANCH, etc.)
            cleaned_model = [
                l for l in model
                if not (l.startswith("ROOT") or l.startswith("ENDROOT") or
                        l.startswith("BRANCH") or l.startswith("ENDBRANCH") or
                        l.strip() == "END")
            ]
            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt") as tmp:
                temp_pdbqt_path = tmp.name

            # Convert/fix ligand model
            fixed_ligand_lines = convert_ligand_model_to_pdb(cleaned_model, temp_pdbqt_path)

            # Write combined protein + fixed ligand
            output_filename = f"{output_base}_model{model_index}.pdb"
            with open(output_filename, 'w') as out:
                out.writelines(protein_lines)
                out.write("\n")
                out.writelines(fixed_ligand_lines)
                out.write("END\n")

            writing_model = False
        elif writing_model:
            model.append(line)

    print(f"Wrote {model_index} combined complex PDB files as {output_base}_model#.pdb.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python combine_and_fix_models.py ligand_out.pdbqt protein.pdb [output_basename]")
        sys.exit(1)

    pdbqt_file = sys.argv[1]
    protein_file = sys.argv[2]
    output_base = sys.argv[3] if len(sys.argv) > 3 else "complex"

    split_and_combine_models(pdbqt_file, protein_file, output_base)

