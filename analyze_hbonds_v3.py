import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import hydrogenbonds
import numpy as np
import matplotlib.pyplot as plt

def analyze_hbonds(topology_file, trajectory_file, prefix, residue_range=None, debug=False, frame_step=1):
    # Load the universe with PSF (topology) and trajectory (coordinates)
    u = mda.Universe(topology_file, trajectory_file)
    
    # Define hydrogen bond analysis with specified criteria
    h = hydrogenbonds.HydrogenBondAnalysis(
        u, 'protein', 'protein', d_h_cutoff=3.5, d_a_cutoff=3.5, d_h_a_angle_cutoff=120.0
    )
    
    # Run analysis with a limited range if in debug mode
    start_frame = 0
    stop_frame = 100 if debug else None
    h.run(start=start_frame, stop=stop_frame, step=frame_step)

    # Prepare output files to save data incrementally
    hbond_detail_file = open(f'{prefix}_hbonds_with_details.txt', 'w')
    hbond_detail_file.write("Donor_resid Donor_resname Donor_name Hydrogen_resid Hydrogen_resname Hydrogen_name "
                            "Acceptor_resid Acceptor_resname Acceptor_name Frame Distance Angle\n")

    hbond_count_file = open(f'{prefix}_hbond_counts_per_frame.txt', 'w')
    hbond_count_file.write("Total hydrogen bonds per frame\n")

    # Prepare atom mapping dictionary
    atom_ids = {atom.id: atom for atom in u.atoms}

    # Prepare variable for summary statistics
    frame_counts = []
    
    # Process hydrogen bonds per frame in chunks
    for ts in u.trajectory[start_frame:stop_frame:frame_step]:
        frame_hbonds = [bond for bond in h.results.hbonds if bond[0] == ts.frame]
        frame_count = 0

        for bond in frame_hbonds:
            # Access donor, hydrogen, and acceptor atoms using atom IDs
            donor_atom = atom_ids.get(bond[1])
            hydrogen_atom = atom_ids.get(bond[2])
            acceptor_atom = atom_ids.get(bond[3])

            # Check that all atoms are valid
            if donor_atom is None or hydrogen_atom is None or acceptor_atom is None:
                continue

            # Filter based on residue range if specified
            if residue_range and not (
                residue_range[0] <= donor_atom.resid <= residue_range[1] or
                residue_range[0] <= acceptor_atom.resid <= residue_range[1]):
                continue

            # Confirm that hydrogen is bonded to nitrogen and acceptor is oxygen
            if donor_atom.name.startswith("N") and hydrogen_atom.name.startswith("H") and acceptor_atom.name.startswith("O"):
                # Write bond details incrementally
                hbond_detail_file.write(
                    f"{donor_atom.resid} {donor_atom.resname} {donor_atom.name} "
                    f"{hydrogen_atom.resid} {hydrogen_atom.resname} {hydrogen_atom.name} "
                    f"{acceptor_atom.resid} {acceptor_atom.resname} {acceptor_atom.name} "
                    f"{ts.frame} {bond[4]:.2f} {bond[5]:.2f}\n"
                )
                frame_count += 1

        # Record frame count and write to file incrementally
        frame_counts.append(frame_count)
        hbond_count_file.write(f"{frame_count}\n")

    # Close the files after writing all data
    hbond_detail_file.close()
    hbond_count_file.close()

    # Plot hydrogen bonds per frame
    plt.plot(range(len(frame_counts)), frame_counts, marker='o', linestyle='-')
    plt.xlabel('Frame')
    plt.ylabel('Number of Hydrogen Bonds')
    plt.title('Hydrogen Bonds per Frame')
    plt.savefig(f'{prefix}_hbonds_per_frame_plot.png')
    plt.show()

    # Histogram of hydrogen bonds per frame
    plt.hist(frame_counts, bins='scott', edgecolor='black')
    plt.xlabel('Number of Hydrogen Bonds')
    plt.ylabel('Frequency')
    plt.title('Histogram of Hydrogen Bonds per Frame')
    plt.savefig(f'{prefix}_hbonds_histogram.png')
    plt.show()

    # Save histogram data and calculate summary statistics
    hist, bin_edges = np.histogram(frame_counts, bins='scott')
    hist_data = np.column_stack((bin_edges[:-1], hist))
    np.savetxt(f'{prefix}_histogram_data.txt', hist_data, fmt='%f', header="Bin Start, Frequency")

    mean_count = np.mean(frame_counts)
    std_count = np.std(frame_counts)
    min_count = np.min(frame_counts)
    max_count = np.max(frame_counts)

    # Write summary statistics to file
    with open(f'{prefix}_summary_statistics.txt', 'w') as f:
        f.write(f"Mean hydrogen bonds per frame: {mean_count:.2f}\n")
        f.write(f"Standard deviation: {std_count:.2f}\n")
        f.write(f"Minimum hydrogen bonds in a frame: {min_count}\n")
        f.write(f"Maximum hydrogen bonds in a frame: {max_count}\n")

def main():
    parser = argparse.ArgumentParser(description='Analyze hydrogen bonds in a trajectory.')
    parser.add_argument('topology_file', help='Path to the topology file.')
    parser.add_argument('trajectory_file', help='Path to the trajectory file.')
    parser.add_argument('prefix', help='Prefix for the output files.')
    parser.add_argument('--residue_range', nargs=2, type=int, default=None,
                        help='Residue range to filter hydrogen bonds (start end).')
    parser.add_argument('--debug', action='store_true', help='Limit analysis to the first 100 frames for debugging.')
    parser.add_argument('--frame_step', type=int, default=1, help='Step size for frames to analyze.')
    args = parser.parse_args()

    analyze_hbonds(args.topology_file, args.trajectory_file, args.prefix, args.residue_range, args.debug, args.frame_step)

if __name__ == '__main__':
    main()

