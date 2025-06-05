"""Command line entry point for the analysis package."""

import argparse
from pathlib import Path

from .clustering import extract_top_clusters, run_hdbscan_clustering
from .pipeline import Pipeline


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

    pipeline_p = subparsers.add_parser("pipeline", help="Run analysis pipeline")
    pipeline_p.add_argument("trajectory", type=Path)
    pipeline_p.add_argument("topology", type=Path)
    pipeline_p.add_argument("prefix", type=str)
    pipeline_p.add_argument("--config", type=Path, default=None)
    pipeline_p.add_argument("--n_clusters", type=int, default=10)
    pipeline_p.add_argument("--full", action="store_true", help="Run full workflow including preprocessing")
    pipeline_p.add_argument("--n_chunks", type=int, default=4)
    pipeline_p.add_argument("--atom_range", type=str, default="1-100")
    pipeline_p.add_argument("--length", type=int, default=1000)
    pipeline_p.add_argument("--merge_selection", type=str, default="all")
    pipeline_p.add_argument("--pca_selection", type=str, default="all")
    pipeline_p.add_argument("--corr_range", type=str, default=None)
    pipeline_p.add_argument("--bins", type=int, default=None)

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
    elif args.command == "pipeline":
        pipeline = Pipeline(args.config)
        if args.full:
            pipeline.full_analysis(
                args.trajectory,
                args.topology,
                args.prefix,
                args.n_chunks,
                args.atom_range,
                args.length,
                n_clusters=args.n_clusters,
                merge_selection=args.merge_selection,
                pca_selection=args.pca_selection,
                corr_range=args.corr_range,
                bins=args.bins,
            )
        else:
            pipeline.cluster_trajectory(
                args.trajectory,
                args.topology,
                args.prefix,
                args.n_clusters,
            )


if __name__ == "__main__":  # pragma: no cover
    main()
