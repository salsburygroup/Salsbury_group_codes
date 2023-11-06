import mdtraj


class Reader:
    def __init__(self, trajectory_path):
        self.trajectory = trajectory_path

    def load(self):
        raise NotImplementedError


class DCD(Reader):
    def __init__(self, trajectory_path, topology_path):
        self.topology = topology_path
        super().__init__(trajectory_path)

    def load(self):
        trajectory = mdtraj.load(self.trajectory, top=self.topology)
        return trajectory


class BigDCD(DCD):
    def __init__(self, trajectory_path, topology_path, chunk_size):
        self.chunk_size = chunk_size
        super().__init__(trajectory_path, topology_path)

    def load(self):
        chunk_generator = mdtraj.iterload(self.trajectory, top=self.topology, chunk=self.chunk_size)
        return chunk_generator


class PDB(Reader):
    def load(self):
        trajectory = mdtraj.load(self.trajectory)
        return trajectory