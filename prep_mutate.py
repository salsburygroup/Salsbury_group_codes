import numpy as np
from simtk.openmm import Vec3
from simtk.openmm.app import PDBFile
from simtk.unit import nanometer, nanometers, Quantity
import pdbfixer
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process PDB file for ACEMD3.')
parser.add_argument('pdbfile', metavar='pdbfile', type=str, help='Input PDB file name')
parser.add_argument('pdbid', metavar='pdbid', type=str, help='Input PDB code')
parser.add_argument('chain', metavar='chain', type=str, help='Chain identifier')
parser.add_argument('resid', metavar='resid', type=int, help='Residue number')
parser.add_argument('mutation', metavar='mutation', type=str, help='Three-letter code of the amino acid to mutate to')
args = parser.parse_args()

# Initialize PDBFixer
fixer = pdbfixer.PDBFixer(filename=args.pdbfile)

# Identify the residue to mutate
mutations = []
for residue in fixer.topology.residues():
    if residue.chain.id == args.chain and residue.id == str(args.resid):
        old_residue_name = residue.name
        new_residue_name = args.mutation
        mutations.append(f"{old_residue_name}-{residue.id}-{new_residue_name}")

# Apply the mutation
fixer.applyMutations(mutations, args.chain)

# Add missing atoms and assign bond orders
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)

# Centering the protein at (0,0,0)
positions = np.array([pos.value_in_unit(nanometers) for pos in fixer.positions])
mean_coords = np.mean(positions, axis=0)
centered_positions = positions - mean_coords
centered_positions = Quantity(centered_positions, nanometers)

# Write the fixed and centered PDB 
out_filename = f"{args.pdbid}_resid{args.resid}_mut{args.mutation}_fixed.pdb"
PDBFile.writeFile(fixer.topology, centered_positions, open(out_filename, 'w'))

