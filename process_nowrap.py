"""Legacy wrapper preserved for backward compatibility."""

from __future__ import annotations

import argparse
from pathlib import Path

from analysis import Pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MD analysis pipeline without wrapping")
    parser.add_argument("trajectory", type=Path)
    parser.add_argument("topology", type=Path)
    parser.add_argument("prefix", type=str)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--n_clusters", type=int, default=10)
    args = parser.parse_args()

    Pipeline(args.config).cluster_trajectory(args.trajectory, args.topology, args.prefix, args.n_clusters)


if __name__ == "__main__":
    main()
