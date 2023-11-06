import mdtraj
import pandas
import numpy

class Processor:
    def __init__(self, trajectory, atom_selection):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        self.sel = atom_selection

    def process(self):
        raise NotImplementedError


class AtomIndexer(Processor):
    # From Oliver Schillinger's mdtraj tools https://github.com/schilli/Tools
    def process(self):
        self.atom_indices = self.trajectory.topology.select(self.sel)
        return self.atom_indices


class Stripper(Processor):
    def process(self):
        raise NotImplementedError


class Strider(Processor):
    def process(self):
        raise NotImplementedError


class Aligner(Processor):
    def __init__(self, trajectory, atom_selection, reference=None):
        self.reference = reference
        super().__init__(trajectory, atom_selection)

    def process(self):
        indices = AtomIndexer(self.trajectory, self.sel).process()
        if self.reference:
            aligned_trajectory = self.trajectory.superpose(reference=self.reference, atom_indices=indices)
        else:
            aligned_trajectory = self.trajectory.superpose(self.trajectory, atom_indices=indices)
        return aligned_trajectory


class Wrapper(Processor):
    def process(self):
        raise NotImplementedError


class Unwrapper(Processor):
    def process(self):
        raise NotImplementedError


class Repair:
    @staticmethod
    def atom_order(right_pdb_file, wrong_pdb_file, wrong_trajectory_file, out_trajectory_file):
        right_topology = mdtraj.load_topology(right_pdb_file)
        wrong_topology = mdtraj.load_topology(wrong_pdb_file)
        right_table, right_bonds = right_topology.to_dataframe()
        wrong_table, wrong_bonds = wrong_topology.to_dataframe()
        repair_list = []
        for index, row in right_table.iterrows():
            repair_list.append(wrong_table[
                            (wrong_table['name'] == row['name'])
                            & (wrong_table['resSeq'] == row['resSeq'])
                            & (wrong_table['resName'] == row['resName'])
                            & (wrong_table['chainID'] == row['chainID'])].index[0])
        with mdtraj.formats.DCDTrajectoryFile(out_trajectory_file, 'w') as right_trajectory_file:
            for frame in mdtraj.iterload(wrong_trajectory_file, top=wrong_pdb_file, chunk=1):
                correct_frame_xyz = frame.xyz[0][repair_list]
                right_trajectory_file.write(correct_frame_xyz*10)  # because nanometers
