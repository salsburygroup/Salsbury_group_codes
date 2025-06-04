"""Pipeline utilities for running analysis steps."""

from __future__ import annotations

import logging
from pathlib import Path

from .config import load_config
from .clustering import extract_top_clusters, run_hdbscan_clustering


class Pipeline:
    """Coordinate execution of analysis steps based on a configuration."""

    def __init__(self, config_path: str | Path | None = None) -> None:
        self.cfg = load_config(config_path)
        self.analysis_dir = Path(self.cfg["analysis_dir"])
        # Paths to external tools are kept for backward compatibility but are
        # no longer used directly by the refactored pipeline.
        self.script_path = Path(self.cfg.get("script_path", "."))
        self.conda_path = Path(self.cfg.get("conda_path", "."))

    def mkdir(self, path: str | Path) -> None:
        Path(path).mkdir(parents=True, exist_ok=True)

    def setup_directories(self) -> None:
        """Create the standard analysis directory layout."""
        for sub in ["SCRATCH", "RMSD", "TRAJ", "RMSF", "CORR", "PCA", "START", "HDBSCAN", "AH", "HBONDS"]:
            self.mkdir(self.analysis_dir / sub)

    def cluster_trajectory(self, traj: Path, top: Path, prefix: str, n_clusters: int) -> None:
        """Run HDBSCAN and extract the most populated clusters."""
        self.setup_directories()
        prefix_path = self.analysis_dir / "HDBSCAN" / prefix
        labels = run_hdbscan_clustering(traj, top, prefix_path)
        extract_top_clusters(
            traj,
            top,
            prefix_path,
            n_clusters,
            prefix_path.with_name(prefix_path.name + "_HDBSCAN_frame_clusters.txt"),
            output_dir=self.analysis_dir / "HDBSCAN",
        )
