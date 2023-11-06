#! /usr/bin/env python
from Analysis import Distance, Saver, TrajectoryReader, TrajectoryProcessor
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# Initialize parser.
parser = argparse.ArgumentParser(
    description='Plot rmsd and rmsf for all simulations', add_help=False)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default=os.getcwd().split('/')[-1])

inputs.add_argument('-n',
                    action='store',
                    dest='number',
                    help='The number of runs',
                    type=int,
                    default=5)

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='name CA')

inputs.add_argument('-tm',
                    action='store',
                    dest='timestep',
                    help='Simulation timestep (ps)',
                    type=int,
                    default=4)

inputs.add_argument('-fq',
                    action='store',
                    dest='dcdfreq',
                    help='dcd frequency',
                    type=int,
                    default='2500')

inputs.add_argument('-dpi',
                    action='store',
                    dest='dpi',
                    help='Dots per inch (resolution of the picture)',
                    type=int,
                    default=200)

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix',
                    type=str,
                    default='CA')

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir=UserInput.out_dir
fname = './' + 'deviation_' + out_dir
if not os.path.exists(fname):
    os.mkdir(fname)

# Xlabel
if UserInput.sel == 'name CA':
    xlabel = 'residue number'
elif UserInput.sel == 'protein': 
    xlabel = 'atom number'
else:
    xlabel = UserInput.sel

# Plot rmsf
fig2=plt.figure(2,figsize=(16,4))
plt.plot(np.loadtxt('/deac/phy/salsburyGrp/wud18/md/TM/TM_thrombin/rmsf_averaged_CA/rmsf_averaged_CA.dat')*10, label='thrombin', alpha=0.9, linewidth=1)
plt.plot(np.loadtxt('/deac/phy/salsburyGrp/wud18/md/TM/TM56_Na0.125/rmsf_averaged_CA/rmsf_averaged_CA.dat')*10, label='TM56-thrombin', alpha=0.9, linewidth=1)
plt.plot(np.loadtxt('/deac/phy/salsburyGrp/wud18/md/TM/TM456_Na0.125/rmsf_averaged_CA/rmsf_averaged_CA.dat')*10, label='TM456-thrombin', alpha=0.9, linewidth=1)

plt.xticks(fontsize=14, rotation=0)
plt.yticks(fontsize=14, rotation=0)
plt.legend(loc='upper right', fontsize=14)
plt.xlabel(xlabel, fontsize=16)
plt.ylabel('RMSF (Å)', fontsize=16)
plt.title('RMSF of ' + UserInput.title, fontsize=20)
fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.png'), pad_inches=0.03, bbox_inches='tight', dpi=UserInput.dpi)
fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.pdf'), bbox_inches='tight')
