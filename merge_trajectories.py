import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align

def concatenate_trajectories(files, topology, output, chunk_size, reference_structure, selection):
    ref = mda.Universe(topology, reference_structure)
    ref_atoms = ref.select_atoms(selection)

    with mda.Writer(f"{output}_merged.xtc", ref_atoms.n_atoms) as W:
        for file in files:
            u = mda.Universe(topology, file)
            u_atoms = u.select_atoms(selection)
            alignment = align.AlignTraj(u, ref, select=selection)
            alignment.run(step=chunk_size)

            for ts in u.trajectory:
                align.alignto(u_atoms, ref_atoms)
                W.write(u_atoms)


def main():
    parser = argparse.ArgumentParser(description="Concatenate and align multiple MD trajectories.")
    parser.add_argument("files", type=str, nargs="+", help="The trajectory files to concatenate.")
    parser.add_argument("--topology", type=str, required=True, help="The topology file.")
    parser.add_argument("--reference_structure", type=str, required=True, help="The reference structure file for alignment.")
    parser.add_argument("--output", type=str, required=True, help="The prefix for the output file.")
    parser.add_argument("--chunk_size", type=int, default=1000, help="The number of frames processed at once.")
    parser.add_argument("--selection", type=str, default="protein", help="The atom selection for alignment.")

    args = parser.parse_args()
    concatenate_trajectories(args.files, args.topology, args.output, args.chunk_size, args.reference_structure, args.selection)

if __name__ == "__main__":
    main()

