"""Command line entry point for the analysis package."""

import argparse
from pathlib import Path

from .clustering import extract_top_clusters, run_hdbscan_clustering
from .pipeline import Pipeline
from .workflow import Workflow


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

    pipeline_p = subparsers.add_parser("pipeline", help="Run simple analysis pipeline")
    pipeline_p.add_argument("trajectory", type=Path)
    pipeline_p.add_argument("topology", type=Path)
    pipeline_p.add_argument("prefix", type=str)
    pipeline_p.add_argument("--config", type=Path, default=None)
    pipeline_p.add_argument("--n_clusters", type=int, default=10)

    full_p = subparsers.add_parser("full", help="Run full workflow")
    full_p.add_argument("prefix", type=str)
    full_p.add_argument("--config", type=Path, default=None)

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
        pipeline.cluster_trajectory(
            args.trajectory,
            args.topology,
            args.prefix,
            args.n_clusters,
        )
    elif args.command == "full":
        wf = Workflow(args.config)
        cmds = wf.build_commands(
            prefix=args.prefix,
            n=4,
            atom_range="1-10",
            path=wf.cfg.get("script_path", "group_python"),
            pdb_prefix="ionized",
            psf_prefix="ionized",
            length=1000,
            nowat_psf_prefix=f"{args.prefix}_autopsf",
            merge_selection="all",
            corr_range=None,
            pca_selection="all",
            bins=None,
            conda_path=wf.cfg.get("conda_path", ""),
            n_clusters=10,
        )
        wf.run_workflow(cmds)


if __name__ == "__main__":  # pragma: no cover
    main()
