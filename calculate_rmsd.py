import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def calculate_rmsd(input_xtc, input_top, input_name):
    # Load the XTC file with the corresponding topology file
    trajectory = md.load_xtc(input_xtc, top=input_top)

    # Get the reference structure (the first frame)
    reference_frame = trajectory[0]

    # Calculate alpha carbon RMSD
    alpha_carbon_indices = trajectory.topology.select('name CA')
    alpha_carbon_rmsd = md.rmsd(trajectory, reference_frame, atom_indices=alpha_carbon_indices)
    alpha_carbon_rmsd_angstrom = alpha_carbon_rmsd * 10  # Convert to angstroms
    np.savetxt(f'{input_name}_alpha_carbon_rmsd.txt', alpha_carbon_rmsd_angstrom)

    # Calculate all-atom RMSD
    all_atom_rmsd = md.rmsd(trajectory, reference_frame)
    all_atom_rmsd_angstrom = all_atom_rmsd * 10  # Convert to angstroms
    np.savetxt(f'{input_name}_all_atom_rmsd.txt', all_atom_rmsd_angstrom)

    # Plot the RMSD values
    time = np.arange(trajectory.n_frames) * trajectory.timestep / 1000  # Time in picoseconds, converted to nanoseconds

    # Plot for Alpha Carbon RMSD
    fig1, ax1 = plt.subplots()
    ax1.plot(time, alpha_carbon_rmsd_angstrom, label=f'{input_name} Alpha Carbon RMSD')
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('RMSD (Å)')
    ax1.legend()
    ax1.set_title(f'{input_name} Alpha Carbon RMSD of Trajectory')
    fig1.savefig(f'{input_name}_alpha_carbon_rmsd.png')

    # Plot for All-Atom RMSD
    fig2, ax2 = plt.subplots()
    ax2.plot(time, all_atom_rmsd_angstrom, label=f'{input_name} All-Atom RMSD')
    ax2.set_xlabel('Time (ns)')
    ax2.set_ylabel('RMSD (Å)')
    ax2.legend()
    ax2.set_title(f'{input_name} All-Atom RMSD of Trajectory')
    fig2.savefig(f'{input_name}_all_atom_rmsd.png')

    # Display both figures
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate the alpha carbon RMSD and all-atom RMSD for a given trajectory file based on the first structure.')
    parser.add_argument('input_xtc', type=str, help='Input XTC file')
    parser.add_argument('input_top', type=str, help='Input topology file (e.g., PDB)')
    parser.add_argument('input_name', type=str, help='Input name to prepend to file names and labels')

    args = parser.parse_args()

    calculate_rmsd(args.input_xtc, args.input_top, args.input_name)

