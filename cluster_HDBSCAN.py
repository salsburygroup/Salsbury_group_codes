# Script 1: cluster_md_part1.py

import argparse
import mdtraj as md
import hdbscan
import numpy as np
import pickle
from collections import Counter
import matplotlib.pyplot as plt

# Parse command line arguments
parser = argparse.ArgumentParser(description='Cluster MD simulations with HDBSCAN (Part 1).')
parser.add_argument('trajectory', type=str, help='Path to the trajectory file.')
parser.add_argument('topology', type=str, help='Path to the topology file.')
parser.add_argument('prefix', type=str, help='Prefix for output file names.')
args = parser.parse_args()

# Load the trajectory
traj = md.load(args.trajectory, top=args.topology)

# Compute pairwise RMSD distances
distances = np.empty((traj.n_frames, traj.n_frames))
for i in range(traj.n_frames):
    distances[i] = md.rmsd(traj, traj, i)

# Perform HDBSCAN clustering
clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=2, metric='precomputed')
labels = clusterer.fit_predict(distances)

# Write out cluster assignments for each frame
with open(f'{args.prefix}_HDBSCAN_frame_clusters.txt', 'w') as f:
    for i, label in enumerate(labels):
        f.write(f'Frame {i}: Cluster {label}\n')

# Compute and write out cluster occupancies, sorted by occupancy
occupancies = Counter(labels)
sorted_occupancies = sorted(occupancies.items(), key=lambda x: x[1], reverse=True)
with open(f'{args.prefix}_HDBSCAN_cluster_occupancies.txt', 'w') as f:
    for label, count in sorted_occupancies:
        f.write(f'Cluster {label}: Occupancy {count}\n')


# Plot cluster number versus time with "+" markers
plt.figure(figsize=(10, 6))
plt.plot(labels, '+')
plt.xlabel('Time (frame number)')
plt.ylabel('Cluster number')
plt.title('Cluster number versus time')
plt.savefig(f'{args.prefix}_HDBSCAN_cluster_vs_time.png')
plt.show(block=True)

# Plot cluster occupancies versus cluster number, sorted by occupancy
plt.figure(figsize=(10, 6))
labels_sorted, counts_sorted = zip(*sorted_occupancies)
plt.bar(range(len(labels_sorted)), counts_sorted, tick_label=labels_sorted)
plt.xlabel('Cluster number')
plt.ylabel('Occupancy')
plt.title('Cluster occupancies versus cluster number')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{args.prefix}_HDBSCAN_occupancy_vs_cluster.png')
plt.show(block=True)

# Plot cumulative sum of sorted occupancies
plt.figure(figsize=(10, 6))
cumulative_counts = np.cumsum(counts_sorted)
plt.plot(range(len(labels_sorted)), cumulative_counts, marker='o')
plt.xlabel('Cluster number')
plt.ylabel('Cumulative occupancy')
plt.title('Cumulative sum of sorted occupancies')
plt.tight_layout()
plt.savefig(f'{args.prefix}_HDBSCAN_cumulative_occupancy.png')
plt.show(block=True)

