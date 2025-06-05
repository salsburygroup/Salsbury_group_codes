"""High level workflow orchestrating all analysis steps."""
from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path
from typing import List

from .pipeline import Pipeline


class Workflow(Pipeline):
    """Run the full analysis workflow using external helper scripts."""

    def _run(self, cmd: List[str] | str, *, shell: bool = False) -> None:
        logging.info("Running: %s", cmd if isinstance(cmd, str) else " ".join(cmd))
        subprocess.run(cmd, check=True, shell=shell)

    def build_commands(
        self,
        prefix: str,
        n: int,
        atom_range: str,
        path: str,
        pdb_prefix: str,
        psf_prefix: str,
        length: int,
        nowat_psf_prefix: str,
        merge_selection: str,
        corr_range: str | None,
        pca_selection: str,
        bins: int | None,
        conda_path: str,
        n_clusters: int,
    ) -> List[List[str] | str]:
        cmds: List[List[str] | str] = []
        script_path = Path(path) / "group_python"
        analysis_dir = Path(self.cfg.get("analysis_dir", "ANALYSIS"))
        scratch_dir = analysis_dir / "SCRATCH"
        rmsd_dir = analysis_dir / "RMSD"
        traj_dir = analysis_dir / "TRAJ"
        rmsf_dir = analysis_dir / "RMSF"
        corr_dir = analysis_dir / "CORR"
        pca_dir = analysis_dir / "PCA"
        start_dir = analysis_dir / "START"

        for d in [analysis_dir, scratch_dir, rmsd_dir, traj_dir, rmsf_dir, corr_dir, pca_dir, start_dir]:
            cmds.append(["mkdir", "-p", str(d)])

        first_value = length // 10
        second_value = length // 100
        if second_value == 1:
            second_value = 2
        first_multiplier = 1000
        second_multiplier = 10000
        smaller_number = min(first_value, second_value)
        larger_number = max(first_value, second_value)

        for i in range(1, n + 1):
            for value, mult in zip([first_value, second_value], [first_multiplier, second_multiplier]):
                cmds.append([
                    sys.executable,
                    str(script_path / "process_trajectory_juststride.py"),
                    f"{prefix}_{i}.xtc",
                    f"{pdb_prefix}.pdb",
                    f"{prefix}_{value}_{i}.xtc",
                    f"{prefix}_{value}_{i}.pdb",
                    str(value),
                    str(mult),
                ])

                cmds.append([
                    conda_path + "conda",
                    "run",
                    "-n",
                    "htmd",
                    "python",
                    str(script_path / "wrap.py"),
                    "-t",
                    f"{prefix}_{value}_{i}.xtc",
                    "-s",
                    f"{psf_prefix}.psf",
                    "-o",
                    f"{prefix}_solvated_wrapped_{value}_{i}.xtc",
                ])

                cmds.append([
                    sys.executable,
                    str(script_path / "process_trajectory_large_atomrange.py"),
                    f"{prefix}_solvated_wrapped_{value}_{i}.xtc",
                    f"{prefix}_{value}_{i}.pdb",
                    f"{prefix}_nowat_{value}_{i}.xtc",
                    f"{prefix}_nowat_{value}_{i}.pdb",
                    "1",
                    str(mult),
                    atom_range,
                ])

                cmds.append([
                    sys.executable,
                    str(script_path / "calculate_rmsd.py"),
                    f"{prefix}_nowat_{value}_{i}.xtc",
                    f"{prefix}_nowat_{value}_{i}.pdb",
                    f"{prefix}_{value}_{i}",
                ])

        for value, mult in zip([first_value, second_value], [first_multiplier, second_multiplier]):
            input_xtcs = [f"{prefix}_nowat_{value}_{i}.xtc" for i in range(1, n + 1)]
            cmds.append([
                sys.executable,
                str(script_path / "merge_trajectories.py"),
                *input_xtcs,
                "--topology",
                f"{nowat_psf_prefix}.psf",
                "--reference_structure",
                f"{prefix}_nowat_{value}_3.pdb",
                "--output",
                f"{prefix}_nowat_{value}",
                "--chunk_size",
                str(mult),
                "--selection",
                merge_selection,
            ])

        cmds.append(["mv", "*_rmsd.txt", str(rmsd_dir) + "/",])
        cmds.append(["mv", "*_rmsd.png", str(rmsd_dir) + "/",])
        cmds.append(["mv", f"{prefix}_nowat_{first_value}_merged.xtc", str(traj_dir) + "/",])
        cmds.append(["mv", f"{prefix}_nowat_{second_value}_merged.xtc", str(traj_dir) + "/",])
        cmds.append(["mv", f"{prefix}_nowat_{smaller_number}_1.pdb", str(traj_dir) + "/",])
        for i in range(1, n + 1):
            cmds.append(["mv", f"{prefix}_{i}.xtc", str(start_dir) + "/",])
        cmds.append(["mv", f"{pdb_prefix}.pdb", str(start_dir) + "/",])
        cmds.append(["mv", f"{psf_prefix}.psf", str(start_dir) + "/",])
        cmds.append(["mv", f"{nowat_psf_prefix}.psf", str(start_dir) + "/",])
        cmds.append(["mv", "*.xtc", str(scratch_dir) + "/",])
        cmds.append(["mv", "*.pdb", str(scratch_dir) + "/",])

        cmds.append([
            sys.executable,
            str(script_path / "calculate_rmsf.py"),
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc",
            prefix,
        ])
        cmds.append(["mv", "*rmsf*", str(rmsf_dir) + "/",])

        if corr_range is not None:
            cmds.append([
                sys.executable,
                str(script_path / "calculate_corr.py"),
                f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
                f"{traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc",
                f"{prefix}_{smaller_number}",
                "--range",
                corr_range,
            ])
        else:
            cmds.append([
                sys.executable,
                str(script_path / "calculate_corr.py"),
                f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
                f"{traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc",
                f"{prefix}_{smaller_number}",
                "--align",
            ])
        cmds.append(["mv", "*correlation_matrix*", str(corr_dir) + "/",])

        cmds.append([
            sys.executable,
            str(script_path / "calculate_pca_projections.py"),
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc",
            pca_selection,
            prefix,
        ])
        cmds.append([
            sys.executable,
            str(script_path / "calculate_pca_corner.py"),
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc",
            pca_selection,
            prefix,
            "5",
        ])

        import math
        bin_number = round(1 + math.log2(n * length * 100 / smaller_number))
        if bins is None:
            bins = bin_number
        cmds.append([
            sys.executable,
            str(script_path / "calculate_FES_PCA.py"),
            f"{prefix}_projections.txt",
            prefix,
            "-b",
            str(bins),
        ])
        cmds.append([
            sys.executable,
            str(script_path / "find_minima_find_structures.py"),
            f"{prefix}_bins_{bins}_free_energy.txt",
            f"{prefix}_bin_indices.txt",
            prefix,
        ])
        cmds.append([
            sys.executable,
            str(script_path / "extract_structures_minima.py"),
            f"{prefix}_projection_minima.txt",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
            prefix,
        ])
        cmds.append(["mv", "*projection*", str(pca_dir) + "/",])
        cmds.append(["mv", "*energy*", str(pca_dir) + "/",])
        cmds.append(["mv", "*minima*", str(pca_dir) + "/",])
        cmds.append(["mv", "*bin*", str(pca_dir) + "/",])
        cmds.append(["mv", "*corner*", str(pca_dir) + "/",])

        cmds.append([
            sys.executable,
            str(script_path / "cluster_HDBSCAN.py"),
            f"{traj_dir}/{prefix}_nowat_{larger_number}_merged.xtc",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
            prefix,
        ])
        cmds.append([
            sys.executable,
            str(script_path / "cluster_traj.py"),
            f"{traj_dir}/{prefix}_nowat_{larger_number}_merged.xtc",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
            prefix,
            str(n_clusters),
            f"{prefix}_HDBSCAN_frame_clusters.txt",
        ])
        hdbscan_dir = analysis_dir / "HDBSCAN"
        cmds.append(["mkdir", "-p", str(hdbscan_dir)])
        cmds.append(["mv", "*HDBSCAN*", str(hdbscan_dir) + "/",])
        cmds.append(["mv", f"{prefix}/*", str(hdbscan_dir) + "/",])
        cmds.append(["rmdir", prefix])

        cmds.append([
            sys.executable,
            str(script_path / "cluster_AH.py"),
            f"{traj_dir}/{prefix}_nowat_{larger_number}_merged.xtc",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
            prefix,
        ])
        cmds.append([
            sys.executable,
            str(script_path / "cluster_traj.py"),
            f"{traj_dir}/{prefix}_nowat_{larger_number}_merged.xtc",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb",
            prefix,
            str(n_clusters),
            f"{prefix}_AH_frame_clusters.txt",
        ])
        ah_dir = analysis_dir / "AH"
        cmds.append(["mkdir", "-p", str(ah_dir)])
        cmds.append(["mv", "*AH*", str(ah_dir) + "/",])
        cmds.append(["mv", f"{prefix}/*", str(ah_dir) + "/",])
        cmds.append(["rmdir", prefix])

        cmds.append([
            sys.executable,
            str(script_path / "HBONDS" / "find_hbonds.py"),
            f"{start_dir}/{nowat_psf_prefix}.psf",
            f"{traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc",
            prefix,
        ])
        cmds.append([
            sys.executable,
            str(script_path / "HBONDS" / "filter_hbonds_cutoff.py"),
            prefix,
        ])
        cmds.append([
            sys.executable,
            str(script_path / "HBONDS" / "cluster_hbonds.py"),
            prefix,
        ])
        hbond_dir = analysis_dir / "HBONDS"
        cmds.append(["mkdir", "-p", str(hbond_dir)])
        cmds.append(["mv", "*hbond*", str(hbond_dir) + "/",])
        cmds.append(["mv", "*info*", str(hbond_dir) + "/",])
        cmds.append(["rmdir", prefix])
        return cmds

    def run_workflow(self, commands: List[List[str] | str]) -> None:
        for cmd in commands:
            self._run(cmd, shell=isinstance(cmd, str))
