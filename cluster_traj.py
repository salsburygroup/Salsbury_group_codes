import argparse
from analysis.clustering import extract_top_clusters


def main() -> None:
    parser = argparse.ArgumentParser(description="Cluster MD simulations with HDBSCAN (Part 2).")
    parser.add_argument("trajectory", type=str, help="Path to the trajectory file.")
    parser.add_argument("topology", type=str, help="Path to the topology file.")
    parser.add_argument("prefix", type=str, help="Prefix for output file names.")
    parser.add_argument("n_clusters", type=int, help="Number of clusters to process.")
    parser.add_argument("labels_file", type=str, help="File with the cluster labels for each frame.")
    args = parser.parse_args()

    extract_top_clusters(args.trajectory, args.topology, args.prefix, args.n_clusters, args.labels_file)


if __name__ == "__main__":
    main()

