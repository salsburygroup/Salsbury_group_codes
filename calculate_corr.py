import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_alpha_carbon_correlation_matrix(xtc_file, pdb_file, align_traj, align_range):
    # Load the trajectory using MDAnalysis
    u = mda.Universe(pdb_file, xtc_file)

    alignment = "none"
    # Align the trajectory if align_traj or align_range is True
    if align_traj:
        reference = u.select_atoms('protein')
        align.AlignTraj(u, reference, select='protein', in_memory=True).run()
        alignment = "protein"
    elif align_range:
        start, end = map(int, align_range.split(":"))
        reference = u.select_atoms(f'bynum {start}:{end}')
        align.AlignTraj(u, reference, select=f'bynum {start}:{end}', in_memory=True).run()
        alignment = f'range_{start}_{end}'

    # Select alpha-carbon atoms
    selection = u.select_atoms('name CA')

    # Get the number of alpha-carbon atoms
    num_atoms = len(selection)

    # Initialize an array to store the positions at each frame
    positions = np.zeros((len(u.trajectory), num_atoms, 3))

    # Iterate over each frame in the trajectory to store positions
    for i, ts in enumerate(u.trajectory):
        positions[i] = selection.positions

    # Calculate mean position for each atom
    mean_positions = positions.mean(axis=0)

    # Calculate positional fluctuations
    fluctuations = positions - mean_positions

    # Compute the covariance matrix
    covariance_matrix = np.einsum('ijk,ilk->jl', fluctuations, fluctuations) / len(u.trajectory)

    # Calculate the standard deviations
    std_devs = np.sqrt(np.diag(covariance_matrix))

    # Normalize the covariance matrix to get the correlation matrix
    correlation_matrix = covariance_matrix / np.outer(std_devs, std_devs)
    np.fill_diagonal(correlation_matrix, 1.0)

    return correlation_matrix, alignment

def plot_correlation_matrix(correlation_matrix, filename, alignment):
    sns.heatmap(correlation_matrix, cmap='coolwarm', center=0)
    plt.title('Alpha-Carbon Correlation Matrix')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.savefig(f"{filename}_{alignment}_correlation_matrix.png")
    plt.show()
    plt.close()

def save_correlation_matrix(correlation_matrix, filename, alignment):
    output_file = f"{filename}_{alignment}_correlation_matrix.txt"
    np.savetxt(output_file, correlation_matrix)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate and plot alpha-carbon correlation matrix from XTC trajectory and PDB topology.')
    parser.add_argument('pdb_file', help='Path to PDB topology file')
    parser.add_argument('xtc_file', help='Path to XTC trajectory file')
    parser.add_argument('output_name', help='Name for the saved correlation matrix')
    parser.add_argument('--align', action='store_true', help='Align the trajectory to the protein before analysis')
    parser.add_argument('--range', type=str, help='Specify an atom range (start:end) for alignment instead of the whole protein')

    args = parser.parse_args()

    correlation_matrix, alignment = calculate_alpha_carbon_correlation_matrix(args.xtc_file, args.pdb_file, args.align, args.range)

    plot_correlation_matrix(correlation_matrix, args.output_name, alignment)

    save_correlation_matrix(correlation_matrix, args.output_name, alignment)

