#! /usr/bin/env python
from Analysis import TrajectoryReader, TrajectoryProcessor
import argparse
import mdtraj

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Align a DCD to reference frame',
    add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Reference PDB',
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

inputs.add_argument('-o',
                    action='store',
                    dest='out_name',
                    help='Output file name',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Read trajectory
reference = TrajectoryReader.PDB(trajectory_path=UserInput.structure).load()
trajectory = TrajectoryReader.BigDCD(
    topology_path=UserInput.structure,
    trajectory_path=UserInput.trajectory,
    chunk_size=1,
).load()

# Align and save frame-wise
with mdtraj.formats.DCDTrajectoryFile(UserInput.out_name, 'w') as output_file:
    for chunk in trajectory:
        aligned_chunk = TrajectoryProcessor.Aligner(
            trajectory=chunk,
            atom_selection=UserInput.sel,
            reference=reference,
        ).process()
        output_file.write(aligned_chunk.xyz*10)
