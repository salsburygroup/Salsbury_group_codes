import sys
import numpy as np
import matplotlib.pyplot as plt

# Ensure that an input filename is provided
if len(sys.argv) < 2:
    print("Usage: python script_name.py input_file.txt")
    sys.exit(1)

input_file = sys.argv[1]

# Load the data, skipping the header row
data = np.loadtxt(input_file, skiprows=1)

# Determine the total number of frames
total_frames = len(data)

# Compute histogram using Scott's rule for bin width
hist, bin_edges = np.histogram(data, bins='fd')

# Convert counts to percentages
hist_percent = (hist / total_frames) * 100.0

# Save histogram bin edges and percentage values to a file
out_data = np.column_stack((bin_edges[:-1], hist_percent))
np.savetxt("hist_data.txt", out_data, header="bin_start  percentage", comments='')

# Plot the histogram as percentages
# Note: We use the previously determined bin_edges instead of 'scott' here
weights = np.ones_like(data) * 100.0 / total_frames

plt.figure(figsize=(6,4))
plt.hist(data, bins=bin_edges, weights=weights, edgecolor='black')
plt.xlabel("Hbonds per Frame")
plt.ylabel("Percentage (%)")
plt.title("Histogram of Hbonds per Frame (Percentage)")
plt.tight_layout()
plt.savefig("histogram.png", dpi=300)
plt.close()

