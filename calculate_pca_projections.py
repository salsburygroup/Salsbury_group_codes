import MDAnalysis as mda
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Handle command line arguments
parser = argparse.ArgumentParser(description='Perform PCA on a molecular dynamics trajectory.')
parser.add_argument('pdb_file', help='The PDB file.')
parser.add_argument('trajectory_file', help='The trajectory file.')
parser.add_argument('selection', help='The selection criteria for the atoms.')
parser.add_argument('output_prefix', help='The prefix for the output files.')
args = parser.parse_args()

# Load the universe
u = mda.Universe(args.pdb_file, args.trajectory_file, in_memory=False)

# Select the atoms
atomgroup = u.select_atoms(args.selection, updating=True)

# Create an empty list to store the coordinates
coordinates = []

# Loop over the trajectory and append the atom coordinates to the list
for ts in u.trajectory:
    coordinates.append(atomgroup.positions)

# Convert the list of coordinates to a numpy array
coordinates = np.array(coordinates)

# Perform PCA on the coordinates
pca = PCA(n_components=2)
projection = pca.fit_transform(coordinates.reshape(coordinates.shape[0], -1))  # Reshape is needed because PCA in scikit-learn does not accept 3D input

# Save the projections to a text file
np.savetxt(f'{args.output_prefix}_projections.txt', projection)

# Plot the projections
plt.scatter(projection[:, 0], projection[:, 1])
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA projection')
plt.savefig(f'{args.output_prefix}_pca_projection.png')

