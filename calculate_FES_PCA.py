import numpy as np
import matplotlib.pyplot as plt
import argparse
import math

# Constants
k_b = 0.0019872041  # Boltzmann constant in kcal/(mol*K)
T = 300  # Temperature in K

# Handle command line arguments
parser = argparse.ArgumentParser(description='Calculate 2D free energy surface from PCA projections.')
parser.add_argument('input_file', help='Input file containing the PCA projections.')
parser.add_argument('output_prefix', help='The prefix for the output files.')
parser.add_argument('-b', '--bins', type=int, default=None, help='Number of bins to use. If not specified, Sturges\' rule is used.')
args = parser.parse_args()

# Load the projections
projections = np.loadtxt(args.input_file)

# Determine the number of bins using Sturges' rule or user input
num_bins = args.bins if args.bins is not None else int(round(1 + np.log2(projections.shape[0])))

# Calculate 2D histogram and bin edges
hist, x_edges, y_edges = np.histogram2d(projections[:, 0], projections[:, 1], bins=num_bins, density=True)

# Calculate bin indices for each projection
bin_indices = np.vstack([np.digitize(projections[:, i], bins=[x_edges, y_edges][i])-1 for i in range(2)]).T

# Save bin indices to a file
np.savetxt(f'{args.output_prefix}_bin_indices.txt', bin_indices, fmt='%d')

# Calculate free energy only for bins with non-zero population
free_energy = np.full(hist.shape, np.inf)  # initialize with infinities
mask = hist > 0  # create a mask for non-zero elements
free_energy[mask] = -k_b * T * np.log(hist[mask])
# Shift energy scale to set global minimum at 0
min_fe = np.min(free_energy[mask])
free_energy[mask] -= min_fe
# Prepare filename
output_filename = f"{args.output_prefix}_bins_{num_bins}_free_energy"
output_filename2 = f"{args.output_prefix}_free_energy"

# Save free energy surface to a text file
np.savetxt(f'{output_filename}.txt', free_energy)

free_energy_transposed = np.transpose(free_energy)
# Plot free energy surface with linear color scale
plt.figure(figsize=(6, 5))
plt.imshow(free_energy_transposed, origin='lower', extent=(x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]))
plt.colorbar(label='Free energy (kcal/mol)')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.savefig(f'{output_filename}.png')  # Save the plot

# Create a filled contour plot of the free energy surface
plt.figure(figsize=(6, 5))
plt.contourf((x_edges[1:] + x_edges[:-1]) / 2, (y_edges[1:] + y_edges[:-1]) / 2, free_energy_transposed, levels=15, cmap='viridis')
plt.colorbar(label='Free energy (kcal/mol)')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.savefig(f'{output_filename2}_contour.png')  # Save the contour plot

