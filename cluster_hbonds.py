import argparse
import numpy as np
import hdbscan
from sklearn.preprocessing import StandardScaler
import os

def compute_hbond_vectors(hbonds):
    hbond_ids = list(map(tuple, hbonds[:, [1, 2, 3]].astype(int)))
    unique_hbond_ids = list(set(hbond_ids))
    hbond_vectors = np.zeros((int(hbonds[:, 0].max() + 1), len(unique_hbond_ids)))

    for hbond in hbonds:
        frame = int(hbond[0])
        hbond_id = tuple(hbond[[1, 2, 3]].astype(int))
        hbond_index = unique_hbond_ids.index(hbond_id)
        hbond_vectors[frame, hbond_index] = 1

    return hbond_vectors, unique_hbond_ids

def cluster_hbonds(prefix, hbond_file=None):
    # Determine the hbond file name and output prefix
    if hbond_file is None:
        hbond_file = f'{prefix}_hbonds_filtered.txt'
        output_prefix = f'{prefix}_hbonds_filtered'
    else:
        output_prefix = os.path.splitext(hbond_file)[0]

    # Load the hydrogen bond data
    hbonds = np.loadtxt(hbond_file)

    # Compute hbond vectors
    hbond_vectors, unique_hbond_ids = compute_hbond_vectors(hbonds)

    # Scale the data
    scaler = StandardScaler()
    hbond_vectors_scaled = scaler.fit_transform(hbond_vectors)

    # Perform HDBSCAN clustering
    clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=2)
    cluster_labels = clusterer.fit_predict(hbond_vectors_scaled)

    # Write the cluster labels to a file
    with open(f'{output_prefix}_cluster_labels.txt', 'w') as f:
        for i, label in enumerate(cluster_labels):
            f.write(f'Frame {i}: Cluster {label}\n')

    # Compute cluster occupancies
    unique, counts = np.unique(cluster_labels, return_counts=True)
    cluster_occupancies = np.asarray((unique, counts)).T

    # Write the cluster occupancies to a file
    np.savetxt(f'{output_prefix}_cluster_occupancies.txt', cluster_occupancies, fmt='%d')

    # Write the hydrogen bonds in each cluster to a file
    for cluster in unique:
        if cluster != -1:  # exclude noise points labelled as -1
            indices = [i for i, x in enumerate(cluster_labels) if x == cluster]
            cluster_hbonds = [unique_hbond_ids[i] for i in indices]
            np.savetxt(f'{output_prefix}_cluster_{cluster}.txt', cluster_hbonds, fmt='%s')

def main():
    parser = argparse.ArgumentParser(description='Cluster frames by hydrogen bonds.')
    parser.add_argument('prefix', help='Prefix for the output files.')
    parser.add_argument('--hbond_file', help='Optional input hydrogen bond file.', default=None)
    args = parser.parse_args()

    cluster_hbonds(args.prefix, args.hbond_file)

if __name__ == '__main__':
    main()

