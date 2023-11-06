from Analysis import Distance, TrajectoryReader, Featurizer
import pyemma.coordinates as coor
import scipy.spatial.distance
import mdtraj as md
import numpy as np
import argparse

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot RMSF', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-t',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='all')

inputs.add_argument('-o1',
                    action='store',
                    dest='out_name1',
                    help='Output name for .pdb',
                    type=str,
                    default='Representative.pdb')

inputs.add_argument('-o2',
                    action='store',
                    dest='out_name2',
                    help='Output name for .dcd',
                    type=str,
                    default='within1sigma.dcd')

# Functions
def rmsd(structure, trajectory, ref_frame):
    # Read trajectory
    trajectory_load = TrajectoryReader.DCD(topology_path=structure, trajectory_path=trajectory).load()
    # Calculate RMSD
    rmsd_timeseries = Distance.RMSD(
        trajectory=trajectory_load, atom_selection=UserInput.sel, reference_frame=ref_frame
    ).calculate()
    return rmsd_timeseries

def median_structure(structure, trajectory):
    trajectory_load = TrajectoryReader.DCD(topology_path=structure, trajectory_path=trajectory).load()
    trajectory_2d = Featurizer.XYZ(trajectory=trajectory_load).extract()
    cluster_string = np.arange(len(trajectory_load))
    cluster_coords = trajectory_2d
    mean = cluster_coords.mean(axis=0)
    distance = [scipy.spatial.distance.euclidean(row, mean) for row in cluster_coords]
    rep = cluster_string[np.argmin(distance)]
    return rep
    
def frames_within_sigma(structure, trajectory):
    RMSD=rmsd(structure, trajectory, representativeFrame)
    sigma=((sum(RMSD**2))/len(RMSD))**0.5
    frames=[]
    for i in range(len(RMSD)):
        if RMSD[i]<=sigma:
            frames.append(i)
    return frames

# Parse into useful form
UserInput = parser.parse_args()

# Loading structure and trajectory
topfile = UserInput.structure
topology = md.load(topfile).topology
trajfile = UserInput.trajectory
feat = coor.featurizer(topfile)
inp = coor.source(trajfile, feat)
print('trajectory length = ',inp.trajectory_length(0))
print('number of dimension = ',inp.dimension())

# Selecting the representative frame and frames within one sigma
representativeFrame=median_structure(topfile, trajfile)
frames=frames_within_sigma(topfile, trajfile)
frames_process=[]
print(len(frames))
for i in range(len(frames)):
    frames_process.append([0,frames[i]])

# Saving the representative frame and frames within one sigma
representative_frame = md.load_frame(trajfile,representativeFrame,top = topfile)
representative_frame.save_pdb(UserInput.out_name1)
coor.save_traj(inp, frames_process, UserInput.out_name2)
