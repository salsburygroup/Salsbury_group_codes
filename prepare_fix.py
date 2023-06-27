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
args = parser.parse_args()

# Initialize PDBFixer
fixer = pdbfixer.PDBFixer(filename=args.pdbfile)

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
PDBFile.writeFile(fixer.topology, centered_positions, open(args.pdbid + '_fixed.pdb', 'w'))

