#! /usr/bin/env python
from scipy.spatial import distance
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis
import argparse

print('Make sure to use the absolute path for output directory(-o)')
# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Compute closest mean distance between Na+ binding loop and Na+', add_help=False)

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
                    default='protein and resid 262:274')

inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Title of the plot',
                    type=str,
                    default='Closest mean distance between Na$\mathregular{^+}$ and 220s loop')

inputs.add_argument('-tm',
                    action='store',
                    dest='time',
                    help='time interval(ps)',
                    type=str,
                    default=10)

inputs.add_argument('-nm',
                    action='store',
                    dest='nm',
                    help='Output name for data',
                    type=str,
                    default='dist_Na_closest_mean_262_274')

inputs.add_argument('-th',
                    action='store',
                    dest='th',
                    help='Threshold for sodium binding',
                    type=str,
                    default='10')

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output folder for data',
                    type=str,
                    required=True)

# Parse into useful form
UserInput = parser.parse_args()

u = MDAnalysis.Universe(UserInput.structure, UserInput.trajectory)
sodiumLoop=u.select_atoms(UserInput.sel)
sodiumGrp=u.select_atoms("resname SOD")
T=len(u.trajectory)
d_traj=np.zeros((T, 1))
for i in range(T):
    u.trajectory[i] #Go to current frame
    d_traj[i]=np.min(np.mean(distance.cdist(sodiumLoop.positions,sodiumGrp.positions,'euclidean'),axis=0))
np.savetxt(UserInput.out_dir + '/' + UserInput.nm + '.dat', d_traj, delimiter=' ')
np.savetxt(UserInput.out_dir + '/' + 'binding_rate' + UserInput.th + '.dat', sum(d_traj<float(UserInput.th))/len(d_traj), delimiter=' ')

# Make time axis1
t=0.001*float(UserInput.time)*np.array(range(1,T+1))

# Plot
fig1=plt.figure(1,figsize=(16,4))
plt.plot(t,d_traj, linewidth=0.3)
plt.title(UserInput.title)
plt.xlabel('time (ns)')
plt.ylabel('Distance (Å)')
fig1.savefig(UserInput.out_dir + '/' + UserInput.nm + '.png', pad_inches=0.03, bbox_inches='tight', dpi=200)
fig1.savefig(UserInput.out_dir + '/' + UserInput.nm + '.tiff', pad_inches=0.03, bbox_inches='tight', dpi=600)
fig1.savefig(UserInput.out_dir + '/' + UserInput.nm + '.pdf', bbox_inches='tight')
