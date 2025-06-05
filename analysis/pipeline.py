"""Pipeline utilities for running analysis steps."""

from __future__ import annotations

import logging
from pathlib import Path
import subprocess
import math
from typing import Iterable

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

    def run_commands(self, commands: Iterable[list[str]]) -> None:
        """Execute a sequence of shell commands."""
        for cmd in commands:
            logging.info("Running %s", " ".join(cmd))
            subprocess.run(cmd, check=True)

    def run_script(self, script_name: str, *args: str) -> list[str]:
        """Build a Python command for one of the helper scripts."""
        script = self.script_path / script_name
        return ["python", str(script), *args]

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
        run_hdbscan_clustering(traj, top, prefix_path)
        extract_top_clusters(
            traj,
            top,
            prefix_path,
            n_clusters,
            prefix_path.with_name(prefix_path.name + "_HDBSCAN_frame_clusters.txt"),
            output_dir=self.analysis_dir / "HDBSCAN",
        )

    def full_analysis(
        self,
        traj: Path,
        top: Path,
        prefix: str,
        n_chunks: int,
        atom_range: str,
        length: int,
        n_clusters: int = 10,
        merge_selection: str = "all",
        pca_selection: str = "all",
        corr_range: str | None = None,
        bins: int | None = None,
    ) -> None:
        """Run the complete analysis workflow including trajectory processing."""

        self.setup_directories()

        first_value = length // 10
        second_value = length // 100
        first_mult = 1000
        second_mult = 10000
        smaller = min(first_value, second_value)
        larger = max(first_value, second_value)

        cmds: list[list[str]] = []

        for i in range(1, n_chunks + 1):
            for value, mult in ((first_value, first_mult), (second_value, second_mult)):
                cmds.append(
                    self.run_script(
                        "process_trajectory_juststride.py",
                        f"{prefix}_{i}.xtc",
                        str(top),
                        f"{prefix}_{value}_{i}.xtc",
                        f"{prefix}_{value}_{i}.pdb",
                        str(value),
                        str(mult),
                    )
                )
                cmds.append(
                    [
                        "python",
                        str(self.script_path / "wrap.py"),
                        "-t",
                        f"{prefix}_{value}_{i}.xtc",
                        "-s",
                        f"{prefix}.psf",
                        "-o",
                        f"{prefix}_solvated_wrapped_{value}_{i}.xtc",
                    ]
                )
                cmds.append(
                    self.run_script(
                        "process_trajectory_large_atomrange.py",
                        f"{prefix}_solvated_wrapped_{value}_{i}.xtc",
                        f"{prefix}_{value}_{i}.pdb",
                        f"{prefix}_nowat_{value}_{i}.xtc",
                        f"{prefix}_nowat_{value}_{i}.pdb",
                        "1",
                        str(mult),
                        atom_range,
                    )
                )
                cmds.append(
                    self.run_script(
                        "calculate_rmsd.py",
                        f"{prefix}_nowat_{value}_{i}.xtc",
                        f"{prefix}_nowat_{value}_{i}.pdb",
                        f"{prefix}_{value}_{i}",
                    )
                )

        for value, mult in ((first_value, first_mult), (second_value, second_mult)):
            xtcs = [f"{prefix}_nowat_{value}_{i}.xtc" for i in range(1, n_chunks + 1)]
            cmds.append(
                self.run_script(
                    "merge_trajectories.py",
                    *xtcs,
                    "--topology",
                    f"{prefix}_autopsf.psf",
                    "--reference_structure",
                    f"{prefix}_nowat_{value}_3.pdb",
                    "--output",
                    f"{prefix}_nowat_{value}",
                    "--chunk_size",
                    str(mult),
                    "--selection",
                    merge_selection,
                )
            )

        cmds.append(
            self.run_script(
                "calculate_rmsf.py",
                f"{prefix}_nowat_{smaller}_1.pdb",
                f"{prefix}_nowat_{smaller}_merged.xtc",
                prefix,
            )
        )

        if corr_range:
            cmds.append(
                self.run_script(
                    "calculate_corr.py",
                    f"{prefix}_nowat_{smaller}_1.pdb",
                    f"{prefix}_nowat_{smaller}_merged.xtc",
                    f"{prefix}_{smaller}",
                    "--range",
                    corr_range,
                )
            )
        else:
            cmds.append(
                self.run_script(
                    "calculate_corr.py",
                    f"{prefix}_nowat_{smaller}_1.pdb",
                    f"{prefix}_nowat_{smaller}_merged.xtc",
                    f"{prefix}_{smaller}",
                    "--align",
                )
            )

        cmds.append(
            self.run_script(
                "calculate_pca_projections.py",
                f"{prefix}_nowat_{smaller}_1.pdb",
                f"{prefix}_nowat_{smaller}_merged.xtc",
                pca_selection,
                prefix,
            )
        )

        cmds.append(
            self.run_script(
                "calculate_pca_corner.py",
                f"{prefix}_nowat_{smaller}_1.pdb",
                f"{prefix}_nowat_{smaller}_merged.xtc",
                pca_selection,
                prefix,
                "5",
            )
        )

        bin_number = round(1 + math.log2(n_chunks * length * 100 / smaller))
        if bins is None:
            bins = bin_number

        cmds.append(
            self.run_script(
                "calculate_FES_PCA.py",
                f"{prefix}_projections.txt",
                prefix,
                "-b",
                str(bins),
            )
        )

        cmds.append(
            self.run_script(
                "find_minima_find_structures.py",
                f"{prefix}_bins_{bins}_free_energy.txt",
                f"{prefix}_bin_indices.txt",
                prefix,
            )
        )

        cmds.append(
            self.run_script(
                "extract_structures_minima.py",
                f"{prefix}_projection_minima.txt",
                f"{prefix}_nowat_{smaller}_merged.xtc",
                f"{prefix}_nowat_{smaller}_1.pdb",
                prefix,
            )
        )

        cmds.append(
            self.run_script(
                "cluster_HDBSCAN.py",
                f"{prefix}_nowat_{larger}_merged.xtc",
                f"{prefix}_nowat_{smaller}_1.pdb",
                prefix,
            )
        )

        cmds.append(
            self.run_script(
                "cluster_traj.py",
                f"{prefix}_nowat_{larger}_merged.xtc",
                f"{prefix}_nowat_{smaller}_1.pdb",
                prefix,
                str(n_clusters),
                f"{prefix}_HDBSCAN_frame_clusters.txt",
            )
        )

        cmds.append(
            self.run_script(
                "cluster_AH.py",
                f"{prefix}_nowat_{larger}_merged.xtc",
                f"{prefix}_nowat_{smaller}_1.pdb",
                prefix,
            )
        )

        cmds.append(
            self.run_script(
                "cluster_traj.py",
                f"{prefix}_nowat_{larger}_merged.xtc",
                f"{prefix}_nowat_{smaller}_1.pdb",
                prefix,
                str(n_clusters),
                f"{prefix}_AH_frame_clusters.txt",
            )
        )

        cmds.append(
            self.run_script(
                "HBONDS/find_hbonds.py",
                f"{prefix}_autopsf.psf",
                f"{prefix}_nowat_{smaller}_merged.xtc",
                prefix,
            )
        )

        cmds.append(self.run_script("HBONDS/filter_hbonds_cutoff.py", prefix))
        cmds.append(self.run_script("HBONDS/cluster_hbonds.py", prefix))

        self.run_commands(cmds)
