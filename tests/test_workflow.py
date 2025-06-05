from analysis.workflow import Workflow


def test_build_commands(tmp_path):
    wf = Workflow()
    cmds = wf.build_commands(
        prefix="sim",
        n=1,
        atom_range="1-2",
        path=".",
        pdb_prefix="ionized",
        psf_prefix="ionized",
        length=1000,
        nowat_psf_prefix="sim_autopsf",
        merge_selection="all",
        corr_range=None,
        pca_selection="all",
        bins=None,
        conda_path="/usr/bin/",
        n_clusters=1,
    )
    assert cmds, "No commands generated"
