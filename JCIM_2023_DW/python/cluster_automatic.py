from Analysis import TrajectoryReader, Featurizer, AtomSelection
from Analysis.Cluster import Clusterer, Saver, Scorer
from collections import Counter
from math import ceil
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Run and score clustering', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', action='store', dest='structure', help='Structure file corresponding to trajectory',type=str, required=True)
inputs.add_argument('-t', action='store', dest='trajectory', help='Trajectory', type=str, required=True)
inputs.add_argument('-sel', action='store', dest='sel', help='Atom selection', type=str, default='not element H')
inputs.add_argument('-o', action='store', dest='out_dir', help='Absolute path of output directory', type=str, required=True)
inputs.add_argument('-m', action='store', dest='method', help='Clustering Method', type=str, required=True)
inputs.add_argument('-mc', action='store', dest='min_cluster_size', help='Minimum Cluster Size', type=str, default=5)
inputs.add_argument('-ms', action='store', dest='min_samples', help='Minimum Samples', type=str, default=5)
inputs.add_argument('-title', action='store', dest='title', help='Title of the timeseries plot', type=str, required=True)
inputs.add_argument('-tm', action='store', dest='timestep', help='Timestep between two frames', type=str, required=True)

# Parse into useful form
UserInput = parser.parse_args()

print('Make sure to use the absolute path for output directory(-o)')
class Predictor:
    def __init__(self, method, dcd_path, topology_path, atom_selection, out_dir, title,  **kwargs):
        self.method = method  # This should be enum
        self.dcd_path = dcd_path
        self.topology_path = topology_path
        self.atom_selection = atom_selection
        self.out_dir = out_dir
        self.title = title
        self.path = os.path.dirname(__file__) + '/' + out_dir
        self.kwargs = kwargs
        self.labels = []
        self.trajectory = []
        self.trajectory_2d = []

    def predict(self):
        # Prepare DCD
        self.trajectory = TrajectoryReader.DCD(topology_path=self.topology_path, trajectory_path=self.dcd_path).load()
        self.trajectory = AtomSelection.Slice(trajectory=self.trajectory, atom_selection=self.atom_selection).select()
        self.trajectory_2d = Featurizer.XYZ(self.trajectory).extract()

        # Compute the differences between adjacent frames
        differences = self.trajectory_2d[1:, :] - self.trajectory_2d[:-1, :]

        # Compute the squared differences
        squared_differences = differences ** 2

        # Sum along the axis representing the coordinates (axis=1) to get squared distance for each frame
        squared_distances = np.sum(squared_differences, axis=1)

        # Take the square root to get the Euclidean distance
        core_distances_array = np.sqrt(squared_distances)

        # Append the last value to the end of the array
        core_distances_array = np.append(core_distances_array, core_distances_array[-1])
        
        del self.trajectory
        # Cluster
        if self.method == 'HDBSCANC':
            clusterer = getattr(Clusterer, self.method)(self.trajectory_2d, minimum_membership=int(UserInput.min_cluster_size), min_samples=int(UserInput.min_samples), core_distance=core_distances_array)
        elif self.method == 'HDBSCAN':
            clusterer = getattr(Clusterer, self.method)(self.trajectory_2d, minimum_membership=int(UserInput.min_cluster_size), min_samples=int(UserInput.min_samples), core_distance=None)
        else:
            clusterer = getattr(Clusterer, self.method)(self.trajectory_2d, **self.kwargs)

        clusterer.fit()
        self.labels = clusterer.labels
#        print(clusterer.labels)
#        print(len(clusterer.labels))
        # Save Timeseries
        Saver.TimeSeries(out_name=os.path.join(self.out_dir, 'timeseries.txt'), labels=clusterer.labels).save()
        # Plot time series
        frames = np.arange(clusterer.labels.shape[0])
        fig1=plt.figure(1,figsize=(8,4))
        plt.scatter(frames, clusterer.labels, marker='+')
        plt.xlabel('Frame(' + UserInput.timestep + ')')
        plt.ylabel('Cluster')
        plt.title(self.method + ' ' + self.title)
        fig1.savefig(os.path.join(out_dir, 'plot_timeseries.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
        fig1.savefig(os.path.join(out_dir, 'plot_timeseries.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
        fig1.savefig(os.path.join(out_dir, 'plot_timeseries.pdf'), bbox_inches='tight')

        # Plot histogram
        if self.method == 'HDBSCAN':
            # Save core distance
            np.savetxt(os.path.join(out_dir, 'Core_distance.txt'), core_distances_array)
            # Graph zero
            fig0=plt.figure(0,figsize=(8,4))
            plt.hist(core_distances_array, bins=30, rwidth=0.9, edgecolor='black')
            plt.title('Histogram of the core distance for ' + UserInput.title)
            plt.xlabel('Core distance (Å)')
            plt.ylabel('Frequency')
            fig0.savefig(os.path.join(out_dir, 'Core_distance.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
            fig0.savefig(os.path.join(out_dir, 'Core_distance.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
            fig0.savefig(os.path.join(out_dir, 'Core_distance.pdf'), bbox_inches='tight')

        ## Graph one
        data=pd.Series(np.loadtxt(os.path.join(self.out_dir, 'timeseries.txt')))
        fig2=plt.figure(2,figsize=(8,4))
        data.plot.hist(bins=200, rwidth=0.9)
        plt.title('Histogram of ' + UserInput.method + ' for ' + UserInput.title)
        plt.xlabel('Clusters')
        plt.ylabel('Number of frames')
        fig2.savefig(os.path.join(out_dir, 'Histogram_raw.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
        fig2.savefig(os.path.join(out_dir, 'Histogram_raw.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
        fig2.savefig(os.path.join(out_dir, 'Histogram_raw.pdf'), bbox_inches='tight')

        ## Graph two
        # Rearranging clusters
        if UserInput.method == 'HDBSCAN':
            delete_noise= "sed '/-1/'d " + os.path.join(self.out_dir, 'timeseries.txt') + " > " + os.path.join(self.out_dir, 'timeseries_no_noise.txt')
            os.system(delete_noise)
            data_no_noise=pd.Series(np.loadtxt(os.path.join(self.out_dir, 'timeseries_no_noise.txt')))
        else:
            data_no_noise=pd.Series(np.loadtxt(os.path.join(self.out_dir, 'timeseries.txt')))            
        num_of_cluster=list(Counter(data_no_noise).values())
        total_num_of_cluster_no_noise=len(data_no_noise)
        total_num_of_cluster=len(data)
        cluster_high2low=[]
        sum_of_cluster_we_draw=[]

        # Pick clusters with big number of configurations
        while np.max(num_of_cluster)/total_num_of_cluster_no_noise > 0.05:
            sum_of_cluster_we_draw.append(np.max(num_of_cluster))
            cluster_high2low.append(round(np.max(num_of_cluster)*100/total_num_of_cluster,2))
            num_of_cluster.remove(np.max(num_of_cluster))
            if len(num_of_cluster) == 0:
                break
        rest=(len(data) - sum(sum_of_cluster_we_draw))
        cluster_high2low.append(round(rest*100/total_num_of_cluster,2))

        # X axis
        index=[]
        name_list = []
        for i in range(len(cluster_high2low)-1):
            index = index + [i]
            name_list = name_list + [i]
            index = index + [len(cluster_high2low)-1]
            name_list = name_list + ['the rest']

        # Y axis
        Yindex=[]
        uplimit=np.max(cluster_high2low)
        if uplimit/10 > 5:
            Yname_list = [0]
            for i in range(ceil(uplimit/10)):
                Yindex = Yindex + [i*10]
                Yname_list = Yname_list + ['0.' + str(i+1)]
        else:
            Yname_list = [0, 0.05]
            for i in range(ceil(uplimit/5)):
                Yindex = Yindex + [i*5]
                Yname_list = Yname_list + ['0.' + str((i+2)*5)]
            
        # Plot
        fig3=plt.figure(3,figsize=(8,4))
        plt.title('Histogram of ' + UserInput.method + ' for ' + UserInput.title)
        plt.xticks(index, name_list)
        plt.yticks(Yindex, Yname_list) 
        plt.xlabel('Clusters')
        plt.ylabel('Probability')
        rects=plt.bar(range(len(cluster_high2low)), cluster_high2low)
        for rect in rects:
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2, height, str(height)+'%', ha='center', va='bottom')
        fig3.savefig(os.path.join(out_dir, 'histogram.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
        fig3.savefig(os.path.join(out_dir, 'histogram.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
        fig3.savefig(os.path.join(out_dir, 'histogram.pdf'), bbox_inches='tight')
        # Score
        silhouette_score = Scorer.Silhouette(labels=clusterer.labels, data=self.trajectory_2d).evaluate()
        Saver.Score(out_name=os.path.join(self.out_dir, 'Silhouette.txt'), score=silhouette_score).save()
        if hasattr(clusterer, 'centers'):
            ch_score = Scorer.CalinskiHarabasz(labels=clusterer.labels, centers=clusterer.centers, data=self.trajectory_2d).evaluate()
            Saver.Score(out_name=os.path.join(self.out_dir, 'CalinskiHarabasz.txt'), score=ch_score).save()
            db_score = Scorer.DaviesBouldin(labels=clusterer.labels, centers=clusterer.centers, data=self.trajectory_2d).evaluate()
            Saver.Score(out_name=os.path.join(self.out_dir, 'DaviesBouldin.txt'), score=db_score).save()

        self.trajectory = TrajectoryReader.DCD(topology_path=self.topology_path, trajectory_path=self.dcd_path).load()

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
       
        full_trajectory = TrajectoryReader.DCD(topology_path=self.topology_path, trajectory_path=self.dcd_path).load()

        Saver.PDB(out_name=os.path.join(self.out_dir, 'clusters'),
                  labels=self.labels, trajectory=full_trajectory).save()
    
if __name__ == "__main__":

    clustering = Predictor(method=UserInput.method,
                           dcd_path=UserInput.trajectory,
                           topology_path=UserInput.structure,
                           atom_selection=UserInput.sel,
                           out_dir=UserInput.out_dir,
                           title=UserInput.title)

    out_dir=UserInput.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    clustering.predict()
    clustering.row_format()
    clustering.extract_pdbs()
