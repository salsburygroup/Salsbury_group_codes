import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

# Handle command line arguments
parser = argparse.ArgumentParser(description='Plot 2D free energy surface from a text file.')
parser.add_argument('input_file', help='Input file containing the free energy data.')
args = parser.parse_args()

# Load the free energy data
free_energy = np.loadtxt(args.input_file)

# Determine non-infinite values
mask = free_energy < np.inf

# Define the corner plot style
sns.set_style('white')
sns.set_context('paper')
plt.rcParams["axes.edgecolor"] = "0.15"
plt.rcParams["axes.linewidth"]  = 1.25

# Create a filled contour plot
contour = plt.contourf(free_energy, origin='lower', levels=15, cmap='viridis', extent=(0, free_energy.shape[0], 0, free_energy.shape[1]))

# Customizing the plot to make it look like corner plots
plt.grid(True, which='both', color='white', linestyle='-', linewidth=1.5)
cbar = plt.colorbar(contour, label='Free energy')
cbar.ax.yaxis.set_tick_params(color='0.15', width=1.25)
plt.setp(plt.getp(cbar.ax.axes, 'yticklines'), color='0.15')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Free Energy Surface')

# Determine output filename
output_filename = os.path.splitext(args.input_file)[0] + '.png'

# Save the plot to a file
plt.savefig(output_filename, dpi=600)

