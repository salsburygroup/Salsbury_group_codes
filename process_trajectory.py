import argparse
import mdtraj as md

def process_trajectory(input_xtc, input_top, output_xtc, output_pdb, stride):
    # Load the XTC file with the corresponding topology file
    trajectory = md.load_xtc(input_xtc, top=input_top)

    # Extract protein and nucleic acid atoms
    complex_atoms = trajectory.topology.select('(protein or nucleic)')
    complex_trajectory = trajectory.atom_slice(complex_atoms)

    # Align the structures
    reference_frame = complex_trajectory[0]
    aligned_complex_trajectory = complex_trajectory.superpose(reference_frame)

    # Save frames at the user-given interval
    subsampled_trajectory = aligned_complex_trajectory[::stride]

    # Write the resulting trajectory to the given XTC file
    subsampled_trajectory.save_xtc(output_xtc)

    # Save a PDB file with just protein and nucleic acid atoms
    reference_frame.save_pdb(output_pdb)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process an MD trajectory.')
    parser.add_argument('input_xtc', type=str, help='Input XTC file')
    parser.add_argument('input_top', type=str, help='Input topology file (e.g., PDB)')
    parser.add_argument('output_xtc', type=str, help='Output XTC file')
    parser.add_argument('output_pdb', type=str, help='Output PDB file with protein and nucleic acid atoms')
    parser.add_argument('stride', type=int, help='Interval to save frames')

    args = parser.parse_args()

    process_trajectory(args.input_xtc, args.input_top, args.output_xtc, args.output_pdb, args.stride)

