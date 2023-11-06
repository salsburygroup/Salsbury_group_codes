#!/bin/bash

# Plot rmsd and rmsf
python /home/wud18/python/Plot_deviation.py -n 8 -title 'thrombin-TM456 complex'

## Concatenate all trajectories and do the analysis
# protein_Na aligned
/home/wud18/.local/lib/vmd/plugins/LINUXAMD64/bin/catdcd5.2/catdcd -o protein_Na.dcd 1/protein_Na.dcd 2/protein_Na.dcd 3/protein_Na.dcd 4/protein_Na.dcd 5/protein_Na.dcd 6/protein_Na.dcd 7/protein_Na.dcd 8/protein_Na.dcd

python /home/wud18/python/align.py -s protein_Na_right_angle.pdb -sel protein -t protein_Na.dcd -o protein_Na_stride1_aligned.dcd && rm protein_Na.dcd

# stride 10 protein aligned
/home/wud18/.local/lib/vmd/plugins/LINUXAMD64/bin/catdcd5.2/catdcd -o protein_Na_stride10_aligned.dcd -stride 10 protein_Na_stride1_aligned.dcd

# stride 100 protein aligned
/home/wud18/.local/lib/vmd/plugins/LINUXAMD64/bin/catdcd5.2/catdcd -o protein_Na_stride100_aligned.dcd -stride 10 protein_Na_stride10_aligned.dcd

# CA aligned
/home/wud18/.local/lib/vmd/plugins/LINUXAMD64/bin/catdcd5.2/catdcd -o CA_stride1_aligned.dcd 1/CA_aligned.dcd 2/CA_aligned.dcd 3/CA_aligned.dcd 4/CA_aligned.dcd 5/CA_aligned.dcd 6/CA_aligned.dcd 7/CA_aligned.dcd 8/CA_aligned.dcd

# Plot averaged rmsf
python /home/wud18/python/rmsf.py -s CA.pdb -t CA_stride1_aligned.dcd -o averaged_CA

# Calculate correlation matrix
python /home/wud18/python/correlation_matrix.py -s CA.pdb -t CA_stride1_aligned.dcd -align

# Calculate the closest mean distance between Na+ and the 220s loop
/home/wud18/bash/NaDistance.sh
