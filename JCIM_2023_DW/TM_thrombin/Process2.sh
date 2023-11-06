#!/bin/bash

# Convert .xtc to .dcd
mdconvert -o output.dcd output.xtc

# Wrap the water and ions
vmd ../input/ionized.psf output.dcd -dispdev text -e /home/wud18/vmd/VMDScripts/Analysis/Wrap.tcl -args -atomsel protein_or_nucleic -outfile output_wrapped.dcd > wrap.log

# Align all frames to the protein
python /home/wud18/python/align.py -s ../input/ionized.pdb -t output_wrapped.dcd -sel protein -o output_wrapped_aligned.dcd > align.log && rm output_wrapped.dcd

# Extract protein_Na and alpha carbon
python /home/wud18/python/extractDCD.py -s ../input/ionized.pdb -t output_wrapped_aligned.dcd -sel 'protein or resname SOD' -o protein_Na.dcd

python /home/wud18/python/extractDCD.py -s ../protein_Na.pdb -t protein_Na.dcd -sel 'name CA and protein' -o CA.dcd

# Align CA 
python /home/wud18/python/align.py -s ../CA.pdb -t CA.dcd -o CA_aligned.dcd && rm CA.dcd

# Calculating RMSF and RMSD for each trajectory
python /home/wud18/python/deviation.py -align -s ../CA.pdb -t CA_aligned.dcd -o CA

# Converting .dcd to .xtc to save space
mdconvert -o output_wrapped_aligned.xtc output_wrapped_aligned.dcd && rm output_wrapped_aligned.dcd
