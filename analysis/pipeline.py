"""Pipeline utilities for running analysis steps."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Iterable

from .config import load_config


class Pipeline:
    """Coordinate execution of analysis steps based on a configuration."""

    def __init__(self, config_path: str | Path | None = None) -> None:
        self.cfg = load_config(config_path)
        self.analysis_dir = Path(self.cfg["analysis_dir"])
        self.script_path = Path(self.cfg["script_path"])
        self.conda_path = Path(self.cfg["conda_path"])

    def run_commands(self, commands: Iterable[list[str]]) -> None:
        for cmd in commands:
            logging.info("Running %s", " ".join(cmd))
            subprocess.run(cmd, check=True)

    def mkdir(self, path: str | Path) -> None:
        Path(path).mkdir(parents=True, exist_ok=True)

    def setup_directories(self) -> None:
        """Create the standard analysis directory layout."""
        for sub in ["SCRATCH", "RMSD", "TRAJ", "RMSF", "CORR", "PCA", "START", "HDBSCAN", "AH", "HBONDS"]:
            self.mkdir(self.analysis_dir / sub)

    def run_script(self, script_name: str, *args: str) -> list[str]:
        script = self.script_path / script_name
        return ["python", str(script), *args]

    def cluster_trajectory(self, traj: Path, top: Path, prefix: str, n_clusters: int) -> None:
        """Run HDBSCAN and extract clusters."""
        self.setup_directories()
        cmds = [
            self.run_script("cluster_HDBSCAN.py", str(traj), str(top), prefix),
            self.run_script("cluster_traj.py", str(traj), str(top), prefix, str(n_clusters), f"{prefix}_HDBSCAN_frame_clusters.txt"),
        ]
        self.run_commands(cmds)
