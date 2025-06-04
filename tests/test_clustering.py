import numpy as np
import mdtraj as md
from analysis.clustering import run_hdbscan_clustering, extract_top_clusters

def _make_traj(tmp_path):
    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue('ALA', chain)
    a1 = top.add_atom('N', element=md.element.nitrogen, residue=res)
    a2 = top.add_atom('CA', element=md.element.carbon, residue=res)
    top.add_bond(a1, a2)
    xyz = np.random.rand(5, 2, 3)
    traj = md.Trajectory(xyz, top)
    traj_file = tmp_path / 'traj.xtc'
    top_file = tmp_path / 'top.pdb'
    traj.save_xtc(traj_file)
    traj[0].save(top_file)
    return traj_file, top_file

def test_run_hdbscan(tmp_path):
    traj_file, top_file = _make_traj(tmp_path)
    labels = run_hdbscan_clustering(
        traj_file,
        top_file,
        str(tmp_path / 'out'),
        min_cluster_size=2,
        min_samples=1,
    )
    assert len(labels) == 5
    assert (tmp_path / 'out_HDBSCAN_frame_clusters.txt').exists()


def test_extract_top_clusters(tmp_path):
    traj_file, top_file = _make_traj(tmp_path)
    labels = np.array([0, 0, 1, 1, -1])
    labels_file = tmp_path / 'labels.txt'
    with open(labels_file, 'w') as fh:
        for i, lbl in enumerate(labels):
            fh.write(f"Frame {i}: Cluster {lbl}\n")
    extract_top_clusters(traj_file, top_file, 'prefix', 1, labels_file,
                         output_dir=tmp_path / 'out')
    assert (tmp_path / 'out/prefix_cluster_0_trajectory.xtc').exists()
    assert (tmp_path / 'out/prefix_cluster_0_median.pdb').exists()
