import mdtraj


class Selector:
    def __init__(self, trajectory, atom_selection):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory
        self.sel = atom_selection

    def select(self):
        raise NotImplementedError


class Indices(Selector):
    def select(self):
        indices = self.trajectory.top.select(self.sel)
        return indices


class Slice(Selector):
    def select(self):
        indices = self.trajectory.top.select(self.sel)
        sub_trajectory = self.trajectory.atom_slice(atom_indices=indices, inplace=False)
        return sub_trajectory
