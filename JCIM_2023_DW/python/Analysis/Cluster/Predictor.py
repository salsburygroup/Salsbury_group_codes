from . import Clusterer
from .. import AtomSelection
from .. import TrajectoryReader
from .. import Featurizer
from . import Plotter
from . import Saver
from . import Scorer
import os
import glob


class Predictor:
    def __init__(self, method, dcd_path, topology_path, atom_selection, out_dir, **kwargs):
        self.method = method  # This should be enum
        self.dcd_path = dcd_path
        self.topology_path = topology_path
        self.atom_selection = atom_selection
        self.out_dir = out_dir
        self.kwargs = kwargs
        self.labels = []
        self.trajectory = []
        self.trajectory_2d = []

    def predict(self):
        # Prepare DCD
        self.trajectory = TrajectoryReader.DCD(trajectory_path=self.dcd_path, topology_path=self.topology_path).load()
        selection = AtomSelection.Slice(trajectory=self.trajectory, atom_selection=self.atom_selection)
        self.trajectory_2d = Featurizer.XYZ(trajectory=selection).extract()
        # Cluster
        clusterer = getattr(Clusterer, self.method)(self.trajectory_2d, **self.kwargs)
        clusterer.fit()
        self.labels = clusterer.labels
        # Save Timeseries
        Saver.TimeSeries(out_name=os.path.join(self.out_dir, 'timeseries.txt'), labels=clusterer.labels).save()
        # Plot time series
        Plotter.TimeSeries(out_name=os.path.join(self.out_dir, 'timeseries.png'), labels=clusterer.labels).plot()
        # Score
        silhouette_score = Scorer.Silhouette(labels=clusterer.labels, data=self.trajectory_2d).evaluate()
        Saver.Score(out_name=os.path.join(self.out_dir, 'silhouette.txt'), score=silhouette_score).save()
        if hasattr(clusterer, 'centers'):
            ch_score = Scorer.CalinskiHarabasz(
                labels=clusterer.labels, centers=clusterer.centers, data=self.trajectory_2d
            ).evaluate()
            Saver.Score(out_name=os.path.join(self.out_dir, 'CalinskiHarabasz.txt'), score=ch_score).save()
            db_score = Scorer.DaviesBouldin(
                labels=clusterer.labels, centers=clusterer.centers, data=self.trajectory_2d
            ).evaluate()
            Saver.Score(out_name=os.path.join(self.out_dir, 'DaviesBouldin.txt'), score=db_score).save()

    def row_format(self):  # Maybe format row
        try:
            self.labels[0]
        except IndexError:
            print('Method predict must be executed first')

        Saver.ClusterFrames(out_name=os.path.join(self.out_dir, 'row_format.txt'), labels=self.labels).save()

    def extract_pdbs(self):
        try:
            self.labels[0]
        except IndexError:
            print('Method predict must be executed first')

        Saver.PDB(
            out_name=os.path.join(self.out_dir, 'clusters'),
            labels=self.labels, trajectory=self.trajectory, atom_selection=self.atom_selection
        ).save()

    def visualize(self):
        for folder in glob.glob(os.path.join(self.out_dir, 'clusters', 'cluster*')):
            Saver.Shadows(
                out_name=folder, middle=os.path.join(folder, 'rep.pdb'), shadow=os.path.join(folder, 'all.pdb')
            ).save()


if __name__ == "__main__":
    import argparse

    # Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
    parser = argparse.ArgumentParser(description='Run and score clustering', add_help=False)

    # List all possible user input
    inputs = parser.add_argument_group('Input arguments')
    inputs.add_argument('-h', '--help', action='help')
    inputs.add_argument('-top', action='store', dest='structure', help='Structure file corresponding to trajectory', type=str, required=True)
    inputs.add_argument('-traj', action='store', dest='trajectory', help='Trajectory', type=str, required=True)
    inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection', type=str, default='not element H')
    inputs.add_argument('-o', action='store', dest='out_dir', help='Output directory', type=str, required=True)
    inputs.add_argument('-method', action='store', dest='method', help='Clustering Method', type=str, required=True)

    # Parse into useful form
    UserInput = parser.parse_args()

    clustering = Predictor(
        method=UserInput.method,
        dcd_path=UserInput.trajectory,
        topology_path=UserInput.structure,
        atom_selection=UserInput.sel,
        out_dir=UserInput.out_dir
    )

    clustering.predict()
    clustering.row_format()
    clustering.extract_pdbs()
    clustering.visualize()
