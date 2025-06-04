"""Command line entry point for the analysis package."""

import argparse
from pathlib import Path

from .clustering import extract_top_clusters, run_hdbscan_clustering


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Utilities for clustering trajectories using HDBSCAN"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    hdbscan_p = subparsers.add_parser("hdbscan", help="Run clustering")
    hdbscan_p.add_argument("trajectory", type=Path)
    hdbscan_p.add_argument("topology", type=Path)
    hdbscan_p.add_argument("prefix", type=str)

    extract_p = subparsers.add_parser("extract", help="Extract top clusters")
    extract_p.add_argument("trajectory", type=Path)
    extract_p.add_argument("topology", type=Path)
    extract_p.add_argument("prefix", type=str)
    extract_p.add_argument("n_clusters", type=int)
    extract_p.add_argument("labels_file", type=Path)

    args = parser.parse_args()
    if args.command == "hdbscan":
        run_hdbscan_clustering(args.trajectory, args.topology, args.prefix)
    elif args.command == "extract":
        extract_top_clusters(
            args.trajectory,
            args.topology,
            args.prefix,
            args.n_clusters,
            args.labels_file,
        )


if __name__ == "__main__":  # pragma: no cover
    main()
