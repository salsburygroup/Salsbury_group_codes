import argparse
import numpy as np
import MDAnalysis as mda
from moleculekit.molecule import Molecule

parser = argparse.ArgumentParser(description='Wrap a large trajectory in chunks', add_help=False)
inputs = parser.add_argument_group('Input arguments')

inputs.add_argument('-p', '--psf',
                    action='store',
                    dest='psf',
                    help='PSF file corresponding to trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-t', '--trajectory',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-sel', '--selection',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='protein')

inputs.add_argument('-o', '--prefix',
                    action='store',
                    dest='prefix',
                    help='Prefix for the output file',
                    type=str,
                    required=True)

inputs.add_argument('-chunk', '--chunk_size',
                    action='store',
                    dest='chunk_size',
                    help='Chunk size for reading trajectory',
                    type=int,
                    default=1000)

UserInput = parser.parse_args()

# Generate output file name based on prefix
out_file = UserInput.prefix + '_wrapped.xtc'

# Initialize MDAnalysis Universe and MoleculeKit Molecule
u = mda.Universe(UserInput.psf, UserInput.trajectory)

# Get the box dimensions from the first frame of the structure file
default_box = u.dimensions

# Create a writer for the output trajectory
with mda.Writer(out_file, u.atoms.n_atoms) as W:
    # Iterate over trajectory in chunks
    for start in range(0, len(u.trajectory), UserInput.chunk_size):
        # Create a FrameIterator for the current chunk
        end = min(start + UserInput.chunk_size, len(u.trajectory))
        chunk_frames = u.trajectory[start:end]

        # Convert chunk to MoleculeKit Molecule
        chunk = Molecule(UserInput.psf)
        chunk_coords = []
        chunk_box = []
        for ts in chunk_frames:
            chunk_coords.append(ts.positions)
            if ts.dimensions is not None:
                chunk_box.append(ts.dimensions[:3])  # Use only the first three values (box lengths)
            else:
                chunk_box.append(default_box[:3])  # Use default box lengths if no box in the frame

        # Debug print: Check the number of frames and boxes
        print(f"Number of frames: {len(chunk_coords)}, Number of boxes: {len(chunk_box)}")

        chunk.coords = np.array(chunk_coords)
        chunk.box = np.array(chunk_box)

        # Debug print: Check the shapes of coords and box
        print(f"Shape of chunk.coords: {chunk.coords.shape}, Shape of chunk.box: {chunk.box.shape}")

        # Perform the wrapping
        chunk.wrap(UserInput.sel)

        # Write wrapped chunk to output trajectory
        for ts in chunk.iterFrames():
            W.write(ts)

