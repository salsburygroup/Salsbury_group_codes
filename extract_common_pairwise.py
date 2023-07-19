import MDAnalysis as mda
import numpy as np
import argparse

def merge_common_atoms(prefix1, prefix2, xtcname1, xtcname2, stepsize):
    # Load the PDB and XTC files
    u1 = mda.Universe(f'{prefix1}_nowat_1_1.pdb', f'{xtcname1}_nowat_1_merged.xtc')
    u2 = mda.Universe(f'{prefix2}_nowat_1_1.pdb', f'{xtcname2}_nowat_1_merged.xtc')

    # Get the atoms from each structure
    atoms1 = [(atom.resid, atom.resname, atom.name) for atom in u1.atoms]
    atoms2 = [(atom.resid, atom.resname, atom.name) for atom in u2.atoms]

    # Find common atoms
    common_atoms = [atom for atom in atoms1 if atom in atoms2]

    # Find indices of the common atoms
    indices_u1 = [i for i, atom in enumerate(atoms1) if atom in common_atoms]
    indices_u2 = [i for i, atom in enumerate(atoms2) if atom in common_atoms]

    # Write out the common atoms to a new PDB file
    u1.atoms[indices_u1].write(f'{prefix1}_{prefix2}_common.pdb')

    # Create a new XTC file with the common atoms
    with mda.Writer(f'{prefix1}_{prefix2}_common_{stepsize}.xtc', len(indices_u1)) as W:
        for ts in u1.trajectory:
            W.write(u1.atoms[indices_u1])
        for ts in u2.trajectory:
            W.write(u2.atoms[indices_u2])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge common atoms from two structures.')
    parser.add_argument('prefix1', type=str, help='Prefix for the first structure files')
    parser.add_argument('prefix2', type=str, help='Prefix for the second structure files')
    parser.add_argument('xtcname1', type=str, help='Name for the first XTC file')
    parser.add_argument('xtcname2', type=str, help='Name for the second XTC file')
    parser.add_argument('stepsize', type=int, help='Step size for the trajectory')

    args = parser.parse_args()

    merge_common_atoms(args.prefix1, args.prefix2, args.xtcname1, args.xtcname2, args.stepsize)

