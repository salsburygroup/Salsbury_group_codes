#! /usr/bin/env python
from Analysis import Distance, Saver, TrajectoryReader, TrajectoryProcessor
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
    description='Calculate, save and plot RMSF' + 'Calculate RMSD time series. Centers structures at' + ' origin before beginning but does not otherwise align.',
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

inputs.add_argument('-ref',
                    action='store',
                    dest='ref_frame',
                    help='Reference Structure',
                    type=int,
                    default=0)

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

inputs.add_argument('-tm',
                    action='store',
                    dest='timestep',
                    help='Simulation timestep (ps)',
                    type=int,
                    default='4')

inputs.add_argument('-fq',
                    action='store',
                    dest='dcdfreq',
                    help='dcd frequency',
                    type=int,
                    default='2500')

inputs.add_argument('-align',
                    action='store_true',
                    help='Align to atom selection before calculating?',)

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default=os.getcwd().split('/')[-1])

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output directory',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir=UserInput.out_dir
fname = './' + 'deviation_' + out_dir
if not os.path.exists(fname):
    os.mkdir(fname)

# Read trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()

if UserInput.align:
    trajectory = TrajectoryProcessor.Aligner(trajectory=trajectory, atom_selection=UserInput.sel).process()

# Calculate RMSD
rmsd_timeseries = Distance.RMSD(
    trajectory=trajectory, atom_selection=UserInput.sel, reference_frame=UserInput.ref_frame
).calculate()

# Save
Saver.Array(array=rmsd_timeseries, out_name=os.path.join('deviation_' + out_dir, 'rmsd_' + out_dir  + '.dat')).save()

# Plot
a=np.arange(trajectory.n_frames)*UserInput.dcdfreq*UserInput.timestep/1000000
fig1=plt.figure(1,figsize=(16,4))
plt.plot(a,10*rmsd_timeseries, alpha=0.9, linewidth=0.8)
plt.xlabel('time (ns)')
plt.ylabel('RMSD (Å)')
plt.title(UserInput.title + ' RMSD timeseries')
fig1.savefig(os.path.join('deviation_'+out_dir,'rmsd_'+out_dir+'.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
#fig1.savefig(os.path.join('deviation_'+out_dir,'rmsd_'+out_dir+'.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
#fig1.savefig(os.path.join('deviation_'+out_dir,'rmsd_'+out_dir+'.pdf'), bbox_inches='tight')

# Calculate RMSF and the mean and std of RMSF
rmsf_array = Distance.RMSF(trajectory=trajectory,
                           atom_selection=UserInput.sel
                           ).calculate()
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
fig2=plt.figure(2,figsize=(16,4))
plt.plot(10*rmsf_array, alpha=0.9, linewidth=0.8)
plt.xlabel(xlabel)
plt.ylabel('RMSF (Å)')
plt.title('RMSF of ' + UserInput.title)
fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.png'), pad_inches=0.03, bbox_inches='tight', dpi=200)
#fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.tiff'), pad_inches=0.03, bbox_inches='tight', dpi=600)
#fig2.savefig(os.path.join('deviation_'+out_dir,'rmsf_'+out_dir+'.pdf'), bbox_inches='tight')

# Save
Saver.Array(array=rmsf_array,
            out_name=os.path.join('deviation_' + out_dir, 'rmsf_' + out_dir + '.dat')).save()
Saver.Array(array=[10*rmsf_mean, 10*rmsf_std],
            out_name=os.path.join('deviation_' + out_dir, 'rmsf_' + out_dir + '_mean_std.dat')).save()
