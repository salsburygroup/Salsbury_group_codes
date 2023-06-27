import argparse
import numpy as np
import mdtraj as md
import os

def process_trajectory(input_xtc, input_top, output_xtc, output_pdb, stride, chunk_size, atom_range):
    complex_trajectory_chunks = []

    # Parsing the atom_range input
    start, end = map(int, atom_range.split('-'))
    atom_indices = list(range(start, end + 1))  # end is inclusive

    for chunk in md.iterload(input_xtc, top=input_top, chunk=chunk_size):
        # Extract the specified range of atoms
        complex_atoms = np.array(atom_indices)

        print(f"complex_atoms: {complex_atoms}, size: {complex_atoms.size}")  # Debugging line

        # Only proceed if we have selected atoms
        if complex_atoms.size > 0:  # Explicitly check size
            chunk_complex = chunk.atom_slice(complex_atoms)
            # Append to the complex_trajectory_chunks
            complex_trajectory_chunks.append(chunk_complex)

    # Join chunks to form the complete complex trajectory
    complex_trajectory = md.join(complex_trajectory_chunks)

    # Align the structures
    reference_frame = complex_trajectory[0]
    aligned_complex = complex_trajectory.superpose(reference_frame, atom_indices=complex_atoms)

    # Write the output XTC and PDB files
    aligned_complex[::stride].save_xtc(output_xtc)
    
    # Write to a temporary PDB file
    temp_pdb = "temp.pdb"
    aligned_complex[0].save_pdb(temp_pdb)

    # Clean up TER lines in the PDB file
    atom_number = 1
    with open(temp_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM'):
                f_out.write("{:6s}{:5d}{}".format('ATOM', atom_number, line[11:]))
                atom_number += 1
            elif line.startswith('TER'):
                f_out.write('TER   \n')  # Added spaces for maintaining column positions
            else:
                f_out.write(line)
    
    # Remove temporary PDB file
    os.remove(temp_pdb)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_xtc')
    parser.add_argument('input_top')
    parser.add_argument('output_xtc')
    parser.add_argument('output_pdb')
    parser.add_argument('stride', type=int)
    parser.add_argument('chunk_size', type=int)
    parser.add_argument('atom_range', help="Range of atoms to select, as 'start-end'")
    args = parser.parse_args()

    process_trajectory(args.input_xtc, args.input_top, args.output_xtc, args.output_pdb, args.stride, args.chunk_size, args.atom_range)

