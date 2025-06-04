"""Helper functions for clustering molecular dynamics trajectories."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable

import hdbscan
import mdtraj as md
import numpy as np


def run_hdbscan_clustering(
    trajectory_file: str | Path,
    topology_file: str | Path,
    prefix: str,
    *,
    min_cluster_size: int = 2,
    min_samples: int = 2,
) -> np.ndarray:
    """Cluster a trajectory with HDBSCAN.

    Parameters
    ----------
    trajectory_file
        Path to the trajectory to cluster.
    topology_file
        Topology for ``trajectory_file``.
    prefix
        Prefix for the output files.
    min_cluster_size
        Minimum cluster size for HDBSCAN.
    min_samples
        Minimum samples for HDBSCAN.

    Returns
    -------
    np.ndarray
        Cluster labels for each frame in ``trajectory_file``.
    """

    traj = md.load(str(trajectory_file), top=str(topology_file))
    distances = np.array(
        [md.rmsd(traj, traj, i) for i in range(traj.n_frames)], dtype=float
    )

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric="precomputed",
    )
    labels = clusterer.fit_predict(distances)

    with open(f"{prefix}_HDBSCAN_frame_clusters.txt", "w") as fh:
        for i, label in enumerate(labels):
            fh.write(f"Frame {i}: Cluster {label}\n")

    occupancies = Counter(labels)
    sorted_occ = sorted(occupancies.items(), key=lambda x: x[1], reverse=True)
    with open(f"{prefix}_HDBSCAN_cluster_occupancies.txt", "w") as fh:
        for label, count in sorted_occ:
            fh.write(f"Cluster {label}: Occupancy {count}\n")

    return labels


def extract_top_clusters(
    trajectory_file: str | Path,
    topology_file: str | Path,
    prefix: str,
    n_clusters: int,
    labels_file: str | Path,
    output_dir: str | Path | None = None,
) -> None:
    """Save trajectories and representative structures for the most populated clusters."""

    out_dir = Path(output_dir or prefix)
    out_dir.mkdir(parents=True, exist_ok=True)

    traj = md.load(str(trajectory_file), top=str(topology_file))

    labels: list[int] = []
    with open(labels_file, "r") as fh:
        for line in fh:
            _, cluster = line.strip().split(":")
            cluster_num = "".join(c for c in cluster.strip() if c.isdigit() or c == "-")
            labels.append(int(cluster_num))
    labels = np.array(labels)

    occupancies = Counter(labels)
    sorted_occ = sorted(occupancies.items(), key=lambda x: x[1], reverse=True)
    top_clusters = [lbl for lbl, _ in sorted_occ[:n_clusters]]

    for cluster_idx in top_clusters:
        cluster_traj = traj[labels == cluster_idx]
        cluster_traj.save(out_dir / f"{prefix}_cluster_{cluster_idx}_trajectory.xtc")

        rmsds = np.array([md.rmsd(cluster_traj, cluster_traj, j) for j in range(cluster_traj.n_frames)])
        median_frame = np.argmin(rmsds.mean(axis=0))
        median_structure = cluster_traj[median_frame]
        median_structure.save(out_dir / f"{prefix}_cluster_{cluster_idx}_median.pdb")

