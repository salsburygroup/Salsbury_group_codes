import numpy as np
import argparse
from scipy.ndimage import gaussian_filter
from scipy.ndimage.filters import minimum_filter
import matplotlib.pyplot as plt

# Handle command line arguments
parser = argparse.ArgumentParser(description='Identify local minima in a 2D free energy surface.')
parser.add_argument('energy_file', help='Input file containing the 2D free energy surface.')
parser.add_argument('indices_file', help='Input file containing the bin indices.')
parser.add_argument('output_prefix', help='The prefix for the output files.')
args = parser.parse_args()

# Load the energy surface and bin indices
free_energy = np.loadtxt(args.energy_file)
bin_indices = np.loadtxt(args.indices_file, dtype=int)


# Generate bin centers filename and load bin centers
bin_centers_file = args.energy_file.replace('_free_energy.txt', '_centers_free_energy.txt')
bin_centers = np.loadtxt(bin_centers_file)

# Check if bin_indices have the right shape, if not, reshape
if bin_indices.ndim == 1:
    bin_indices = bin_indices.reshape(-1, 2)

# Find the maximum finite value in the free energy surface
max_fe = np.nanmax(free_energy[np.isfinite(free_energy)])

# Replace infinities with a value slightly higher than the maximum finite value
free_energy[np.isinf(free_energy)] = max_fe * 1.01  # Increase by 1%

# Smooth the energy surface with a Gaussian filter to reduce the number of spurious minima
smooth_energy = gaussian_filter(free_energy, sigma=1)

# Save the smoothed energy surface to a text file
np.savetxt(f'{args.output_prefix}_smooth_energy.txt', smooth_energy)

#transpose due to imshow
#smooth_energy_transposed = np.transpose(smooth_energy)

# Plot the smoothed energy surface and save it as a PNG image
plt.imshow(smooth_energy, origin='lower')
plt.colorbar(label='Smoothed energy')
plt.savefig(f'{args.output_prefix}_smooth_energy.png')

# Define an 8-connected neighborhood (this can be easily adjusted to other connectivity if needed)
neighborhood = np.ones((3, 3))

# Apply the minimum filter to the data, then subtract the result from the data
# This leaves only the local minima in the image
min_filtered = minimum_filter(smooth_energy, footprint=neighborhood)
diff = smooth_energy - min_filtered
peaks = np.where(diff == 0)

# Extract the coordinates of the local minima
minima = np.column_stack(peaks) # Assign unique labels to the minima
labels = np.arange(len(minima)) + 1

# Output the big coordinates and labels of the minima
minima_bins = np.column_stack([minima, labels])
np.savetxt(f'{args.output_prefix}_minima_bins.txt', minima_bins, fmt='%d')

# Identify which bins belong to local minima
minima_mask = np.full(bin_indices.shape[0], -1, dtype=int)
minima_coords = minima

projection_minima = []
for i, bin_coord in enumerate(bin_indices):
    for j, minima_coord in enumerate(minima_coords):
        if np.array_equal(bin_coord, minima_coord):
            minima_mask[i] = labels[j]
            projection_minima.append((i + 1, labels[j], minima_coord[0], minima_coord[1]))  # i + 1 because of 1-based indexing and added minima_coord for coordinates
            break

# Output the labels of the local minima that each bin belongs to
projection_minima = np.array(projection_minima)
np.savetxt(f'{args.output_prefix}_projection_minima.txt', projection_minima, fmt='%d')

# Calculate the percentage of total projections for each minima
unique, counts = np.unique(projection_minima[:, 1], return_counts=True)
percentages = counts / bin_indices.shape[0] * 100

# Reshape unique array to 2D
unique = unique.reshape(-1, 1)

# Fetch corresponding minima_coords and bin_centers for unique labels
minima_coords_unique = [minima_coords[np.where(labels == u)[0][0]] for u in unique]
bin_centers_unique = [bin_centers[np.ravel_multi_index(minima_coords[np.where(labels == u)[0][0]], free_energy.shape)] for u in unique]

# Convert lists to 2D arrays
minima_coords_unique = np.array(minima_coords_unique)
bin_centers_unique = np.array(bin_centers_unique)

# Reshape percentages to 2D
percentages = percentages.reshape(-1, 1)

# Join unique, minima_coords_unique, bin_centers_unique, and percentages
minima_percentages_with_coords_centers = np.hstack((unique, minima_coords_unique, bin_centers_unique, percentages))

np.savetxt(f'{args.output_prefix}_minima_percentages.txt', minima_percentages_with_coords_centers, fmt=['%d', '%d', '%d', '%.2f', '%.2f', '%.2f'])

