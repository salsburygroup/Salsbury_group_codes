import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from analysis.clustering import run_hdbscan_clustering


def main() -> None:
    parser = argparse.ArgumentParser(description="Cluster MD simulations with HDBSCAN.")
    parser.add_argument("trajectory", type=str, help="Path to the trajectory file.")
    parser.add_argument("topology", type=str, help="Path to the topology file.")
    parser.add_argument("prefix", type=str, help="Prefix for output file names.")
    args = parser.parse_args()

    labels = run_hdbscan_clustering(args.trajectory, args.topology, args.prefix)

    occupancies = Counter(labels)
    sorted_occ = sorted(occupancies.items(), key=lambda x: x[1], reverse=True)

    plt.figure(figsize=(10, 6))
    plt.plot(labels, "+")
    plt.xlabel("Time (frame number)")
    plt.ylabel("Cluster number")
    plt.title("Cluster number versus time")
    plt.savefig(f"{args.prefix}_HDBSCAN_cluster_vs_time.png")

    plt.figure(figsize=(10, 6))
    labels_sorted, counts_sorted = zip(*sorted_occ)
    plt.bar(range(len(labels_sorted)), counts_sorted, tick_label=labels_sorted)
    plt.xlabel("Cluster number")
    plt.ylabel("Occupancy")
    plt.title("Cluster occupancies versus cluster number")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{args.prefix}_HDBSCAN_occupancy_vs_cluster.png")

    plt.figure(figsize=(10, 6))
    cumulative_counts = np.cumsum(counts_sorted)
    plt.plot(range(len(labels_sorted)), cumulative_counts, marker='o')
    plt.xlabel("Cluster number")
    plt.ylabel("Cumulative occupancy")
    plt.title("Cumulative sum of sorted occupancies")
    plt.tight_layout()
    plt.savefig(f"{args.prefix}_HDBSCAN_cumulative_occupancy.png")


if __name__ == "__main__":
    main()

