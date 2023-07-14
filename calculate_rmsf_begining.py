import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import sys

def calculate_rmsf(structure_file, xtc_file, output_name, selection='all', selection_name=None):
    # If no selection_name is given, use the selection as the name
    if selection_name is None:
        selection_name = selection

    # Load the Universe
    u = mda.Universe(structure_file, xtc_file)

    # Select atoms
    atoms = u.select_atoms(selection)

# Calculate the RMSF
    rmsf = np.zeros((len(atoms),))
    average_coordinates = np.zeros((len(atoms), 3))

    n_frames = min(3000, len(u.trajectory))  # Ensure not to exceed the total number of frames
    start_frame = len(u.trajectory) - n_frames

# Calculate the average coordinates over the last 1000 frames
    for ts in u.trajectory[start_frame:]:
      average_coordinates += atoms.positions
    average_coordinates /= n_frames

# Reset the trajectory iterator to the start of the last 1000 frames
    u.trajectory[start_frame]

# Calculate the rmsf over the last 1000 frames
    for ts in u.trajectory[start_frame:]:
      rmsf += np.sum((atoms.positions - average_coordinates)**2, axis=1)
    rmsf = np.sqrt(rmsf / n_frames)

# Save the RMSF to a file
    np.savetxt(f'{output_name}_{selection_name}_rmsf.txt', rmsf)


    # Plot the RMSF and save the plot
    plt.figure()
    plt.plot(rmsf)
    plt.xlabel('Atom Number')
    plt.ylabel('RMSF (Å) ')
    plt.title(f'RMSF Plot ({selection_name})')
    plt.savefig(f'{output_name}_{selection_name}_rmsf_plot.png')
    plt.close()

    # Write B-factors to PDB file
    atoms_temp = atoms.atoms.tempfactors
    atoms.atoms.tempfactors = rmsf
    with mda.Writer(f'{output_name}_{selection_name}_rmsf.pdb', atoms.n_atoms) as W:
        W.write(atoms)
    atoms.atoms.tempfactors = atoms_temp  # reset the original B-factors

    # Return the RMSF for subplot
    return rmsf

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rmsf.py <structure_file> <xtc_file> <output_name>")
        sys.exit(1)
    structure_file = sys.argv[1]
    xtc_file = sys.argv[2]
    output_name = sys.argv[3]

    # Calculate RMSF for all atoms
    rmsf_all = calculate_rmsf(structure_file, xtc_file, output_name, 'all')

    # Calculate RMSF for alpha-carbon atoms
    rmsf_ca = calculate_rmsf(structure_file, xtc_file, output_name, 'name CA', 'CA')

    # Display the RMSFs in subplots
    fig, axs = plt.subplots(2)
    axs[0].plot(rmsf_all)
    axs[0].set_title('RMSF Plot (all atoms)')
    axs[1].plot(rmsf_ca)
    axs[1].set_title('RMSF Plot (CA atoms)')
    for ax in axs:
        ax.set(xlabel='Atom', ylabel='RMSF')
    plt.tight_layout()
    #plt.show()

