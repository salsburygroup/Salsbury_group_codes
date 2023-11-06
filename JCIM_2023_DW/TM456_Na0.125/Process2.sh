#!/bin/bash

# Extract first 1us simulations
#/home/wud18/.local/lib/vmd/plugins/LINUXAMD64/bin/catdcd5.2/catdcd -o output.dcd -first 1 -last 100000 -xtc output.xtc

# Wrap the water and ions
#vmd ../input/ionized.psf output.dcd -dispdev text -e /home/wud18/vmd/VMDScripts/Analysis/Wrap.tcl -args -atomsel protein_or_nucleic -outfile output_wrapped.dcd > wrap.log

#vmd ../input/ionized.psf output_wrapped.dcd -dispdev text -e /home/wud18/vmd/VMDScripts/Analysis/Wrap.tcl -args -atomsel protein_or_nucleic -outfile output_wrapped_TM.dcd > wrap.log

# Align all frames to the protein
python /home/wud18/python/align.py -s ../input/ionized.pdb -t output_wrapped.dcd -sel protein -o output_wrapped_aligned.dcd > align.log 

python /home/wud18/python/align.py -s ../input/ionized.pdb -t output_wrapped_TM.dcd -sel protein -o output_wrapped_TM_aligned.dcd > align_TM.log 

# Extract protein_Na and alpha carbon
python /home/wud18/python/extractDCD.py -s ../input/ionized.pdb -t output_wrapped_aligned.dcd -sel 'segname AP1 BP1 or resname SOD' -o protein_Na.dcd

python /home/wud18/python/extractDCD.py -s ../protein_Na.pdb -t protein_Na.dcd -sel 'name CA and protein' -o CA.dcd

# Align CA 
python /home/wud18/python/align.py -s ../CA.pdb -t CA.dcd -o CA_aligned.dcd && rm CA.dcd

# Calculating RMSF and RMSD for each trajectory
python /home/wud18/python/deviation.py -align -s ../CA.pdb -t CA_aligned.dcd -o CA

# Converting .dcd to .xtc to save space
mdconvert -o output_wrapped_aligned.xtc output_wrapped_aligned.dcd && rm output_wrapped_aligned.dcd

mdconvert -o output_wrapped.xtc output_wrapped.dcd && rm output_wrapped.dcd

mdconvert -o output_wrapped_TM.xtc output_wrapped_TM.dcd && rm output_wrapped_TM.dcd
