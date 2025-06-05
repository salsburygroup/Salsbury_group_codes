from pathlib import Path
import numpy as np
import mdtraj as md
from analysis.pipeline import Pipeline

from tests.test_clustering import _make_traj


def test_pipeline_cluster(tmp_path):
    traj_file, top_file = _make_traj(tmp_path)
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text(f"analysis_dir: {tmp_path / 'ANALYSIS'}")
    pipeline = Pipeline(cfg_file)
    pipeline.cluster_trajectory(traj_file, top_file, 'prefix', 1)
    files = list((tmp_path / 'ANALYSIS' / 'HDBSCAN').glob('prefix_cluster_*_median.pdb'))
    assert files, "No median structure written"
