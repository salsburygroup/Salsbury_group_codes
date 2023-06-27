import numpy as np
import MDAnalysis as mda
import argparse

# Handle command line arguments
parser = argparse.ArgumentParser(description='Extract frames corresponding to each minima from a trajectory.')
parser.add_argument('projection_minima_file', help='Input file containing projection minima information.')
parser.add_argument('trajectory_file', help='Input trajectory file.')
parser.add_argument('topology_file', help='Input topology file.')
parser.add_argument('output_prefix', help='The prefix for the output files.')
args = parser.parse_args()

# Load projection minima
projection_minima = np.loadtxt(args.projection_minima_file, dtype=int)

# Load the trajectory
u = mda.Universe(args.topology_file, args.trajectory_file)

minima_ids = np.unique(projection_minima[:, 1])
for minima_id in minima_ids:
    # Create a new trajectory writer for each minima
    with mda.Writer(f'{args.output_prefix}_minima_{minima_id}.xtc', u.atoms.n_atoms) as W:
        # Select frames that belong to the current minima
        frames = projection_minima[projection_minima[:, 1] == minima_id, 0] - 1  # 0-indexed

        # Get coordinates for each frame and write them to the new trajectory file
        coords = []
        for frame in frames:
            u.trajectory[frame]  # update universe to the desired frame
            W.write(u)  # write the current TimeStep in the universe
            coords.append(u.atoms.positions)  # save positions for finding median later

        coords = np.array(coords)

        if coords.size > 0:
            # Find the frame that is most central among all frames in the current minima
            median_frame = frames[np.argmin(np.sum(np.square(coords - coords.mean(axis=0)), axis=(1,2)))]

            # Write the most central frame to a separate PDB file
            u.trajectory[median_frame]  # update universe to the median frame
            u.atoms.write(f'{args.output_prefix}_minima_{minima_id}_median.pdb')
        else:
            print(f"No coordinates for minima {minima_id}. Skipping.")

