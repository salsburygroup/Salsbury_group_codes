#!/usr/bin/env python
import numpy as np
import argparse
import mdtraj

# Jiajie Xiao
# Jan 14, 2017

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Generate a new trajectory based on selection',
    add_help=False)

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
                    help='selection',
                    type=str,
                    default='all')

inputs.add_argument('-l',
                    action='store',
                    dest='stride',
                    help='stride to use. It will be ignored when frameFile flas is specified.',
                    type=int,
                    default=1)

inputs.add_argument('-f',
                    action='store',
                    dest='frameFile',
                    help='file with frame of interest',
                    type=str,
                    default='noFile')

inputs.add_argument('-o',
                    action='store',
                    dest='outputfile',
                    help='output dcd',
                    type=str,
                    default='selected.dcd')

# Parse into useful form
UserInput = parser.parse_args()

if UserInput.frameFile == 'noFile':
    traj = mdtraj.load(UserInput.trajectory, top=UserInput.structure, stride=UserInput.stride)
else:
    temp = mdtraj.load(UserInput.trajectory, top=UserInput.structure)
    frame = np.loadtxt(UserInput.frameFile,delimiter=',') 
    traj = temp.slice(frame.astype(int)-1, copy = True)

#print(traj.topology.to_dataframe())
sel = traj.topology.select(UserInput.sel)
extractTraj = traj.atom_slice(sel)

extractTraj.save(UserInput.outputfile)
