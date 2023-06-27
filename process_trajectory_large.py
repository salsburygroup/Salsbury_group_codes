import argparse
import mdtraj as md
import numpy as np

def process_trajectory(input_xtc, input_top, output_xtc, output_pdb, stride, chunk_size):
    # Initialize an empty list for chunks of the complex trajectory
    complex_trajectory_chunks = []

    # Load the XTC file in chunks
    for chunk in md.iterload(input_xtc, top=input_top, chunk=chunk_size):
        # Extract protein and nucleic acid atoms separately
        protein_atoms = chunk.topology.select('protein')
      #  print(f"Protein atoms: {protein_atoms}")
        nucleic_atoms = chunk.topology.select('resn DA DG DC DT')
      #  print(f"Nucleic atoms: {nucleic_atoms}")

        # Combine protein and nucleic atoms
        complex_atoms = np.concatenate((protein_atoms, nucleic_atoms))

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
    aligned_complex[0].save_pdb(output_pdb)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_xtc')
    parser.add_argument('input_top')
    parser.add_argument('output_xtc')
    parser.add_argument('output_pdb')
    parser.add_argument('stride', type=int)
    parser.add_argument('chunk_size', type=int)
    args = parser.parse_args()

    process_trajectory(args.input_xtc, args.input_top, args.output_xtc, args.output_pdb, args.stride, args.chunk_size)

