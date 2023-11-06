#!/usr/bin/env python
from Analysis import AtomSelection, Correlation, Plotter, Saver, TrajectoryReader, TrajectoryProcessor
import mdtraj as md
import numpy as np
import argparse
import os

# Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Calculate, save and plot correlation matrix', add_help=False)

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

inputs.add_argument('-tau',
                    action='store',
                    dest='covariance_tau',
                    default=None,
                    type=int,
                    help='Lag time for constructing a time-lagged correlation matrix',)

inputs.add_argument('-align',
                    action='store_true',
                    help='Align to atom selection before calculating?',)

inputs.add_argument('-ax',
                    action='store',
                    dest='axis_label',
                    default=None,
                    help='Label for axes',)

inputs.add_argument('-directional',
                    action='store',
                    dest='directional',
                    help='directional correlation?yes or no',
                    type=str,
                    default='no')

inputs.add_argument('-th',
                    action='store',
                    dest='th',
                    help='Threshold for correlation coefficients',
                    type=float,
                    default=0.8)

inputs.add_argument('-ig',
                    action='store',
                    dest='ig',
                    help='Ignore adjacent residues',
                    type=int,
                    default=1)

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
                    default='corr')

# Parse into useful form
UserInput = parser.parse_args()

# Make output directory
out_dir=UserInput.out_dir
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Process trajectory
trajectory = TrajectoryReader.DCD(topology_path=UserInput.structure, trajectory_path=UserInput.trajectory).load()
trajectory = AtomSelection.Slice(trajectory=trajectory, atom_selection=UserInput.sel).select()

if UserInput.align:
    trajectory = TrajectoryProcessor.Aligner(trajectory=trajectory, atom_selection=UserInput.sel).process()

# Xlabel
if UserInput.sel == 'name CA':
    xlabel = 'residue number'
elif UserInput.sel == 'protein': 
    xlabel = 'atom number'
else:
    xlabel = UserInput.sel

# Make correlation matrix
if UserInput.covariance_tau:
    correlation_matrix = Correlation.TimeLagged(
        trajectory=trajectory, covariance_tau=UserInput.covariance_tau
    ).calculate()
    title = 'Correlation Matrix with tau = {0}'.format(UserInput.covariance_tau) + ' (' + UserInput.title + ')'
else:
    correlation_matrix = Correlation.Pearson(trajectory=trajectory).calculate()
    title = 'Correlation Matrix (' + UserInput.title + ')'

# Plot HeatMap
if UserInput.axis_label:
    Plotter.UnityPColor(y=correlation_matrix,
                        out_name=os.path.join(out_dir, out_dir),
                        x_label=UserInput.axis_label,
                        y_label=UserInput.axis_label,
                        title=title).plot()
else:
    Plotter.UnityPColor(y=correlation_matrix,
                        out_name=os.path.join(out_dir, out_dir),
                        x_label=xlabel,
                        y_label=xlabel,
                        title=title).plot()

# Save HeatMap
Saver.Array(
    array=correlation_matrix,
    out_name=os.path.join(out_dir, out_dir + '.dat')).save()

topology = md.load(UserInput.structure).topology
residue_list = [residue for residue in topology.residues]
corr = np.genfromtxt(os.path.join(out_dir, out_dir + '.dat'))
num_residues = corr.shape[0]

# Correlation coefficients correspond to residue pairs
with open(os.path.join(out_dir, out_dir+'_table.csv'), 'w') as out:
    out.write('Res_A, Res_B, Type, Corr, abs(Corr)\n')
    for id_a in range(num_residues):
        if UserInput.directional == 'no':
            id_b_start = id_a
        else:
            id_b_start = 0
        for id_b in range(id_b_start, num_residues):
            out.write('{0:s}, {1:s}, {2:s}, {3:f}, {4:f}\n'.format(str(residue_list[id_a]), str(residue_list[id_b]), str(int(np.sign(corr[id_a, id_b]))), corr[id_a, id_b], abs(corr[id_a, id_b])))

# Save non_diagonal values
with open(os.path.join(out_dir, out_dir + '_non_diagonal.dat'), 'w') as out:
    for id_a in range(num_residues):
        id_b_start = id_a + 1
        for id_b in range(id_b_start, num_residues):
            out.write('{0:f}\n'.format(corr[id_a, id_b]))

# Save the average value of non_diagonal values
corr_non_diagonal = np.genfromtxt(os.path.join(out_dir, UserInput.out_dir + '_non_diagonal.dat'))
corr_non_diagonal_mean = corr_non_diagonal.mean()
corr_non_diagonal_square = corr_non_diagonal**2
corr_non_diagonal_square_mean = corr_non_diagonal_square.mean()
corr_non_diagonal_square_mean_sqrt = corr_non_diagonal_square_mean**0.5

Saver.Array(array=[corr_non_diagonal_mean, corr_non_diagonal_square_mean_sqrt],
            out_name=(os.path.join(out_dir, UserInput.out_dir+'_non_diagonal_mean&square_mean_sqrt.dat'))).save()

# Correlation coefficients map to structure
with open(os.path.join(out_dir, out_dir+'_map.dat'), 'w') as out:
    for id_a in range(num_residues):
        id_b_start = id_a
        for id_b in range(id_b_start, num_residues):
            out.write('{0:s}\t{1:s}\t{2:f}\n'.format(str(id_a), str(id_b), corr[id_a, id_b]))

# Correlation coefficients map to structure(corr > UserInput.th)
with open(os.path.join(out_dir, out_dir + '_map_' + str(UserInput.th) + '_ig' + str(UserInput.ig) + '.dat'), 'w') as out:
    for id_a in range(num_residues):
        id_b_start = id_a + UserInput.ig
        for id_b in range(id_b_start, num_residues):
            if abs(corr[id_a, id_b]) > UserInput.th:
                out.write('{0:s}\t{1:s}\t{2:f}\n'.format(str(id_a), str(id_b), corr[id_a, id_b]))
