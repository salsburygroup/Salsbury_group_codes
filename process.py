"""High-level orchestration script for MD analysis."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from analysis import Pipeline


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MD analysis pipeline")
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("topology", type=Path)
    parser.add_argument("prefix", type=str)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--n_clusters", type=int, default=10)
    args = parser.parse_args()

    pipeline = Pipeline(args.config)
    pipeline.cluster_trajectory(args.trajectory, args.topology, args.prefix, args.n_clusters)


if __name__ == "__main__":
    main()
