import argparse
import MDAnalysis as mda

def process_trajectory(input_xtc, input_top, output_xtc, output_pdb, stride, chunk_size):
    # Create a universe from the input XTC and topology files
    universe = mda.Universe(input_top, input_xtc)

    # Write the first frame to the output PDB file
    universe.atoms.write(output_pdb)

    # Write to the output XTC file
    with mda.Writer(output_xtc, universe.atoms.n_atoms) as W:
        for frame in range(0, len(universe.trajectory), stride):
            universe.trajectory[frame]
            print(f"Processing frame: {universe.trajectory.frame}, stride: {stride}")
            W.write(universe)
            chunk_size -= 1
            if chunk_size == 0:
                break

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

