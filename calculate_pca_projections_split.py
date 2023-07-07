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
parser.add_argument('ranges', nargs='+', type=int, help='List of frame number ranges.')
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

# Split the projections based on the ranges
ranges = np.array(args.ranges)
for i in range(0, len(ranges), 2):
    start = ranges[i]
    end = ranges[i+1]
    split_projection = projection[start:end]

    # Save the split projections to a text file
    np.savetxt(f'{args.output_prefix}_projections_{start}_{end}.txt', split_projection)

    # Plot the split projections
    plt.scatter(split_projection[:, 0], split_projection[:, 1])
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title(f'PCA projection for frames {start} to {end}')
    plt.savefig(f'{args.output_prefix}_pca_projection_{start}_{end}.png')
    plt.clf()  # Clear the figure for the next plot

