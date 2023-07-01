import os
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Find the most chemically diverse molecules.')
parser.add_argument('size', type=int, help='The number of input molecules.')
parser.add_argument('number', type=int, help='The number of output molecules.')
args = parser.parse_args()

# Load molecules
molecules = []
file_names = []
for file in os.listdir(f'top_hits{args.size}'):
    if file.startswith('zinc_hitnumber_') and file.endswith('_rdkit.pdb'):
        mol = Chem.MolFromPDBFile(os.path.join(f'top_hits{args.size}', file))
        if mol is not None:
            molecules.append(mol)
            file_names.append(file)

# Calculate fingerprints
fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in molecules]

# Calculate similarity using Tanimoto coefficient
similarity_matrix = []
for i in range(len(fingerprints)):
    similarity_row = []
    for j in range(len(fingerprints)):
        similarity = DataStructs.FingerprintSimilarity(fingerprints[i], fingerprints[j])
        similarity_row.append(similarity)
    similarity_matrix.append(similarity_row)

# Convert similarity to distance
dist_matrix = 1 - np.array(similarity_matrix)

# Find most diverse molecules using MaxMin algorithm
most_diverse_indices = []
most_diverse_indices.append(np.argmax(dist_matrix.sum(axis=1)))
for _ in range(1, args.number):
    remaining_indices = np.setdiff1d(np.arange(len(molecules)), most_diverse_indices)
    min_distances = dist_matrix[most_diverse_indices, :][:, remaining_indices].min(axis=0)
    most_diverse_indices.append(remaining_indices[np.argmax(min_distances)])

# Save most diverse molecules
os.makedirs(f'most_diverse_{args.number}_from{args.size}', exist_ok=True)
for i in most_diverse_indices:
    mol = molecules[i]
    file_name = file_names[i].replace('_rdkit.pdb', '')
    Chem.MolToPDBFile(mol, os.path.join(f'most_diverse_{args.number}_from{args.size}', f'{file_name}.pdb'))
    Chem.MolToMolFile(mol, os.path.join(f'most_diverse_{args.number}_from{args.size}', f'{file_name}.sdf'))

    # Generate SMILES string and write to .smi file
    smiles = Chem.MolToSmiles(mol)
    with open(os.path.join(f'most_diverse_{args.number}_from{args.size}', f'{file_name}.smi'), 'w') as f:
        f.write(smiles)

# Calculate distance matrix of most diverse molecules
most_diverse_distance_matrix = dist_matrix[np.ix_(most_diverse_indices, most_diverse_indices)]

# Save distance matrix to .txt file
np.savetxt(os.path.join(f'most_diverse_{args.number}_from{args.size}', 'distance_matrix.txt'), most_diverse_distance_matrix, fmt='%.6f')

