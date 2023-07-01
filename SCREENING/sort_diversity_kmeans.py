import os
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.cluster import KMeans
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

# Find most diverse molecules using k-means clustering
kmeans = KMeans(n_clusters=args.number)
kmeans.fit(dist_matrix)
cluster_centers = kmeans.cluster_centers_

# Write cluster information to a file
output_dir = f'most_diverse_{args.number}_from{args.size}_kmeans'
os.makedirs(output_dir, exist_ok=True)
with open(os.path.join(output_dir, 'cluster_info.txt'), 'w') as f:
    for i, center in enumerate(cluster_centers):
        indices = [idx for idx, label in enumerate(kmeans.labels_) if label == i]
        center_index = indices[np.argmin(np.linalg.norm(dist_matrix[indices] - center, axis=1))]
        center_molecule = file_names[center_index].replace('_rdkit.pdb', '')
        size = list(kmeans.labels_).count(i)
        f.write(f'Cluster number: {i}, Size: {size}, Center: {center_molecule}\n')

# Save most diverse molecules
most_diverse_indices = [np.argmin(np.linalg.norm(dist_matrix - center, axis=1)) for center in cluster_centers]
for i in most_diverse_indices:
    mol = molecules[i]
    file_name = file_names[i].replace('_rdkit.pdb', '')
    Chem.MolToPDBFile(mol, os.path.join(output_dir, f'{file_name}.pdb'))
    Chem.MolToMolFile(mol, os.path.join(output_dir, f'{file_name}.sdf'))

    # Generate SMILES string and write to .smi file
    smiles = Chem.MolToSmiles(mol)
    with open(os.path.join(output_dir, f'{file_name}.smi'), 'w') as f:
        f.write(smiles)

# Calculate distance matrix of most diverse molecules
most_diverse_distance_matrix = dist_matrix[np.ix_(most_diverse_indices, most_diverse_indices)]

# Save distance matrix to .txt file
np.savetxt(os.path.join(output_dir, 'distance_matrix.txt'), most_diverse_distance_matrix, fmt='%.6f')

