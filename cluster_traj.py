import argparse
import mdtraj as md
import numpy as np
import os
from collections import Counter

# Parse command line arguments
parser = argparse.ArgumentParser(description='Cluster MD simulations with HDBSCAN (Part 2).')
parser.add_argument('trajectory', type=str, help='Path to the trajectory file.')
parser.add_argument('topology', type=str, help='Path to the topology file.')
parser.add_argument('prefix', type=str, help='Prefix for output file names.')
parser.add_argument('n_clusters', type=int, help='Number of clusters to process.')
parser.add_argument('labels_file', type=str, help='File with the cluster labels for each frame.')
args = parser.parse_args()

# Create an output directory
output_dir = args.prefix
os.makedirs(output_dir, exist_ok=True)

# Load the trajectory
traj = md.load(args.trajectory, top=args.topology)

# Load the labels from the text file
labels = []
with open(args.labels_file, 'r') as f:
    for line in f:
        _, cluster = line.strip().split(":")
        cluster_num = ''.join(c for c in cluster.strip() if c.isdigit() or c == '-')  # Strip non-digit characters except '-'
        labels.append(int(cluster_num))

labels = np.array(labels)

# Load the cluster occupancies from the previous script and sort them
occupancies = Counter(labels)
sorted_occupancies = sorted(occupancies.items(), key=lambda x: x[1], reverse=True)

# Get the labels of the first n clusters
n_clusters = args.n_clusters
top_clusters = [label for label, _ in sorted_occupancies[:n_clusters]]

# Write out a trajectory for each of the top clusters
for i in top_clusters:
    cluster_traj = traj[labels == i]
    cluster_traj.save(os.path.join(output_dir, f'{args.prefix}_cluster_{i}_trajectory.xtc'))

# Compute and write out median structure for each of the top clusters
for i in top_clusters:
    cluster_traj = traj[labels == i]
    rmsds = np.empty((cluster_traj.n_frames, cluster_traj.n_frames))
    for j in range(cluster_traj.n_frames):
        rmsds[j] = md.rmsd(cluster_traj, cluster_traj, j)
    median_frame = np.argmin(rmsds.mean(axis=0))
    median_structure = cluster_traj[median_frame]
    median_structure.save(os.path.join(output_dir, f'{args.prefix}_cluster_{i}_median.pdb'))

