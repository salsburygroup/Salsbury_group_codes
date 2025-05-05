import sys
import os

def split_pdbqt_models(input_file, output_base="model"):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    model = []
    model_index = 0
    writing_model = False

    for line in lines:
        if line.startswith("MODEL"):
            writing_model = True
            model = [line]
        elif line.startswith("ENDMDL"):
            model.append(line)
            model_index += 1
            # Remove unwanted PDBQT-specific lines
            cleaned_model = [
                l for l in model
                if not (l.startswith("ROOT") or l.startswith("ENDROOT") or
                        l.startswith("BRANCH") or l.startswith("ENDBRANCH"))
            ]
            output_filename = f"{output_base}_model{model_index}.pdb"
            with open(output_filename, 'w') as out:
                out.writelines(cleaned_model)
            writing_model = False
        elif writing_model:
            model.append(line)

    print(f"Extracted {model_index} model(s) to {output_base}_model#.pdb files.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python split_pdbqt_models.py ligand_out.pdbqt [output_basename]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_base = sys.argv[2] if len(sys.argv) > 2 else "model"
    split_pdbqt_models(input_file, output_base)

