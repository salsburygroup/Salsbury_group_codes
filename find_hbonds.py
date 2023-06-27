import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import hydrogenbonds
import numpy as np

def analyze_hbonds(topology_file, trajectory_file, prefix):
    u = mda.Universe(topology_file, trajectory_file)

    # Define the analysis
    h = hydrogenbonds.HydrogenBondAnalysis(u, 'protein', 'protein', d_h_cutoff=1.2, d_a_cutoff=3.2, d_h_a_angle_cutoff=120.0)

    # Analyze the first 100 frames
    #h.run(start=0, stop=100)
    h.run()
    total_timesteps = len(u.trajectory)
    #total_timesteps = 100
# Access the results
    hbonds = h.results.hbonds
    counts_by_ids = h.count_by_ids()
    counts_by_time = h.count_by_time()

    # Write the results to files with prefix
    np.savetxt(f'{prefix}_hbonds.txt', hbonds, fmt='%s')
    np.savetxt(f'{prefix}_hbond_counts_by_ids.txt', counts_by_ids, fmt='%s')
    np.savetxt(f'{prefix}_hbond_counts_by_time.txt', counts_by_time, fmt='%s')

    # Write out residue name, number, atom name, and atom id for each atom
    atom_info = [[atom.id, atom.resname, atom.resid, atom.name] for atom in u.atoms]
    np.savetxt(f'{prefix}_atom_info.txt', atom_info, fmt='%s')

def main():
    parser = argparse.ArgumentParser(description='Find hydrogen bonds in a trajectory.')
    parser.add_argument('topology_file', help='Path to the topology file.')
    parser.add_argument('trajectory_file', help='Path to the trajectory file.')
    parser.add_argument('prefix', help='Prefix for the output files.')
    args = parser.parse_args()

    analyze_hbonds(args.topology_file, args.trajectory_file, args.prefix)

if __name__ == '__main__':
    main()

