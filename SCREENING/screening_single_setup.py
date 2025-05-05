import sys
import os
import argparse
import numpy as np
import shutil

def compute_box(pdb_file, margin=5.0):
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue
    if not coords:
        raise ValueError("No atomic coordinates found in PDB file.")
    
    coords = np.array(coords)
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    center = (min_coords + max_coords) / 2
    size = (max_coords - min_coords) + 2 * margin
    return center.tolist(), size.tolist()

def generate_scripts(pdb_file, project_name, ligand_name, ligand_file_path, margin):
    pdb_prefix = os.path.splitext(os.path.basename(pdb_file))[0]
    docking_dir = f"/deac/phy/salsburyGrp/autodock/{project_name}/ligand_{ligand_name}"
    ligand_ext = os.path.splitext(ligand_file_path)[1].lower()

    center, size = compute_box(pdb_file, margin)
    center_x, center_y, center_z = [f"{c:.2f}" for c in center]
    size_x, size_y, size_z = [f"{s:.2f}" for s in size]

    os.makedirs("SCREENING", exist_ok=True)

    with open("SCREENING/ligand_prep.sh", "w") as f:
        f.write('#!/bin/bash\n')
        f.write(f'mkdir -p {docking_dir}\n')
        f.write(f'cp {ligand_file_path} {docking_dir}/ligand{ligand_ext}\n')
        f.write(f'cd {docking_dir}\n')
        if ligand_ext == ".mol2":
            f.write(f'/home/salsbufr/babel/bin/babel -i mol2 ligand.mol2 -o pdb ligand.pdb\n')
        else:
            f.write(f'cp ligand{ligand_ext} ligand.pdb\n')  # fallback copy if already PDB
        f.write('/home/luy/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ligand.pdb\n')

    with open("SCREENING/run_vina.sh", "w") as f:
        f.write('#!/bin/bash -l\n')
        f.write('. /etc/profile.d/modules.sh\n')
        f.write('module load apps/vina/1.2.5\n')
        f.write(f'cd {docking_dir}\n')
        f.write(f'vina --receptor {pdb_prefix}.pdbqt --ligand ligand.pdbqt '
                f'--center_x {center_x} --center_y {center_y} --center_z {center_z} '
                f'--size_x {size_x} --size_y {size_y} --size_z {size_z} --cpu 1 > log.txt\n')

    with open(f"docking_{ligand_name}.slurm", "w") as f:
        f.write('#!/bin/bash -l\n')
        f.write('#SBATCH --partition=small\n')
        f.write(f'#SBATCH --job-name="dock_{ligand_name}"\n')
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --ntasks=1\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write('#SBATCH --time=2:00:00\n')
        f.write('#SBATCH --account="salsburyGrp"\n')
        f.write(f'#SBATCH --output="SLURM_OUTS/dock_{ligand_name}_%j.o"\n')
        f.write(f'#SBATCH --error="SLURM_ERRORS/dock_{ligand_name}_%j.e"\n')
        f.write('#SBATCH --mem=1gb\n')
        f.write('source /etc/profile\n')
        f.write(f'mkdir -p {docking_dir}\n')
        f.write(f'cp {pdb_file} {docking_dir}/{pdb_prefix}.pdb\n')
        f.write(f'/home/luy/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py '
                f'-r {docking_dir}/{pdb_prefix}.pdb -o {docking_dir}/{pdb_prefix}.pdbqt\n')
        f.write('bash SCREENING/ligand_prep.sh\n')
        f.write('bash SCREENING/run_vina.sh\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate docking scripts for a single ligand.")
    parser.add_argument("pdb_file", help="Full path to receptor PDB file")
    parser.add_argument("project_name", help="Project name")
    parser.add_argument("ligand_name", help="Ligand label")
    parser.add_argument("ligand_file_path", help="Full path to ligand file (e.g., .pdb, .mol2)")
    parser.add_argument("--margin", type=float, default=5.0, help="Padding in Å added to each box dimension (default: 5.0)")
    args = parser.parse_args()

    generate_scripts(args.pdb_file, args.project_name, args.ligand_name, args.ligand_file_path, args.margin)

