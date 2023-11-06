import mdtraj


class FeatureExtractor:
    def __init__(self, trajectory):
        assert isinstance(trajectory, mdtraj.Trajectory)
        self.trajectory = trajectory

    def extract(self):
        raise NotImplementedError


class XYZ(FeatureExtractor):
    def extract(self):
        xyz_2d = self.trajectory.xyz.reshape(self.trajectory.n_frames, self.trajectory.n_atoms * 3)
        xyz_2d = xyz_2d.astype('float64')
        return xyz_2d


class COM(FeatureExtractor):
    def extract(self):
        com_traj = mdtraj.compute_center_of_mass(self.trajectory)
        return com_traj
