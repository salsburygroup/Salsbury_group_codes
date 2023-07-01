import os
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import hdbscan

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Find the most chemically diverse molecules.')
parser.add_argument('size', type=int, help='The number of input molecules.')
parser.add_argument('num_clusters', type=int, help='Number of the most occupied clusters. -1 for all clusters.')
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

# Find most diverse molecules using HDBSCAN clustering
clusterer = hdbscan.HDBSCAN(min_cluster_size=2)
cluster_labels = clusterer.fit_predict(dist_matrix)

if args.num_clusters == -1:
    # use all clusters
    cluster_centers = [dist_matrix[cluster_labels == i].mean(axis=0) for i in np.unique(cluster_labels) if i != -1]
else:
    # use top n clusters
    cluster_sizes = np.bincount(cluster_labels[cluster_labels!=-1]) # Count number of molecules in each cluster
    top_clusters = np.argsort(cluster_sizes)[-args.num_clusters:]
    cluster_centers = [dist_matrix[cluster_labels == i].mean(axis=0) for i in top_clusters]

most_diverse_indices = [np.argmin(np.linalg.norm(dist_matrix - center, axis=1)) for center in cluster_centers]

# Save most diverse molecules
os.makedirs(f'most_diverse_from{args.size}_hdbscan_{args.num_clusters}', exist_ok=True)
for i in most_diverse_indices:
    mol = molecules[i]
    file_name = file_names[i].replace('_rdkit.pdb', '')
    Chem.MolToPDBFile(mol, os.path.join(f'most_diverse_from{args.size}_hdbscan_{args.num_clusters}', f'{file_name}.pdb'))
    Chem.MolToMolFile(mol, os.path.join(f'most_diverse_from{args.size}_hdbscan_{args.num_clusters}', f'{file_name}.sdf'))

    # Generate SMILES string and write to .smi file
    smiles = Chem.MolToSmiles(mol)
    with open(os.path.join(f'most_diverse_from{args.size}_hdbscan_{args.num_clusters}', f'{file_name}.smi'), 'w') as f:
        f.write(smiles)

# Calculate distance matrix of most diverse molecules
most_diverse_distance_matrix = dist_matrix[np.ix_(most_diverse_indices, most_diverse_indices)]

# Save distance matrix to .txt file
np.savetxt(os.path.join(f'most_diverse_from{args.size}_hdbscan_{args.num_clusters}', 'distance_matrix.txt'), most_diverse_distance_matrix, fmt='%.6f')

# Get number of molecules in each cluster and the center molecule
cluster_sizes = np.bincount(cluster_labels[cluster_labels!=-1]) # Count number of molecules in each cluster
cluster_center_molecules = [file_names[i].replace('_rdkit.pdb', '') for i in most_diverse_indices] # Get center molecules

# Save cluster information to .txt file
with open(os.path.join(f'most_diverse_from{args.size}_hdbscan_{args.num_clusters}', 'cluster_info.txt'), 'w') as f:
    for i, (size, center) in enumerate(zip(cluster_sizes, cluster_center_molecules)):
        f.write(f'Cluster {i+1}: Size = {size}, Center = {center}\n')

