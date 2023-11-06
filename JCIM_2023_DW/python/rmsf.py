#!/usr/bin/env python
from Analysis import Distance, Saver, TrajectoryReader, TrajectoryProcessor
import matplotlib.pyplot as plt
import argparse
import os

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='text-davinci-003 model', add_help=False)

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
                    default='name CA')

inputs.add_argument('-align',
                    action='store',
                    dest='align',
                    help='Align to atom selection before calculating?',
                    type=str,
                    default='name CA')

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default=os.getcwd().split('/')[-1])

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix for text and png',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir=UserInput.out_dir
fname = './rmsf_' + out_dir
if not os.path.exists(fname):
    os.mkdir(fname)

# Load trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = TrajectoryProcessor.Aligner(trajectory=trajectory, atom_selection=UserInput.align).process()

# Calculate RMSF and the mean and std of RMSF
rmsf_array = Distance.RMSF(trajectory=trajectory, atom_selection=UserInput.sel).calculate()
rmsf_mean = rmsf_array.mean()
rmsf_std = rmsf_array.std()

# Xlabel
if UserInput.sel == 'name CA':
    xlabel = 'residue number'
elif UserInput.sel == 'protein': 
    xlabel = 'atom number'
else:
    xlabel = UserInput.sel

# Plot
fig1=plt.figure(1,figsize=(16,4))
plt.plot(10*rmsf_array)
plt.xlabel(xlabel)
plt.ylabel('RMSF (Å)')
plt.title('RMSF of ' + UserInput.title)
fig1.savefig(os.path.join('rmsf_'+out_dir,'rmsf_'+out_dir+'.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
fig1.savefig(os.path.join('rmsf_'+out_dir,'rmsf_'+out_dir+'.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig1.savefig(os.path.join('rmsf_'+out_dir,'rmsf_'+out_dir+'.pdf'), bbox_inches='tight')

# Save
Saver.Array(array=rmsf_array,
            out_name=os.path.join('rmsf_'+out_dir, 'rmsf_'+out_dir+'.dat')).save()
Saver.Array(array=[10*rmsf_mean, 10*rmsf_std],
            out_name=os.path.join('rmsf_'+out_dir, 'rmsf_'+out_dir+'_mean_std.dat')).save()
