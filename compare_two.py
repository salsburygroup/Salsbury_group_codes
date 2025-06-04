import subprocess
import sys
import os

# Getting the prefixes from the command line arguments
prefix1 = sys.argv[1]
prefix2 = sys.argv[2]
script_path = sys.argv[3]

# Define the command templates
cmds = [
    "python {0}/calculate_rmsf_diff.py {1}_all_rmsf.pdb {2}_all_rmsf.pdb {1}_{2}_diff_rmsf.pdb",
    "python {0}/calculate_matrix_difference.py {1}_1_protein_correlation_matrix.txt {2}_1_protein_correlation_matrix.txt {1}_{2}_difference_correlation_matrix",
    "python {0}/extract_common_pairwise.py {1} {2} {1} {2} 1",
    "python {0}/calculate_pca_projections_split.py {1}_{2}_common.pdb {1}_{2}_common_1.xtc all {1}_{2}_common_1 0 39999 40000 79999",
    "python {0}/calculate_FES_PCA.py {1}_{2}_common_1_projections_40000_79999.txt {1}_{2}_{2}",
    "python {0}/calculate_FES_PCA.py {1}_{2}_common_1_projections_0_39999.txt {1}_{2}_{1}",
    "python {0}/find_minima_find_structures.py {1}_{2}_{1}_bins_16_free_energy.txt {1}_{2}_{1}_bin_indices.txt {1}_{2}_{1}",
    "python {0}/find_minima_find_structures.py {1}_{2}_{2}_bins_16_free_energy.txt {1}_{2}_{2}_bin_indices.txt {1}_{2}_{2}",
    "python {0}/extract_structures_minima.py {1}_{2}_{1}_projection_minima.txt {1}_{2}_common_1.xtc {1}_{2}_common.pdb {1}_{2}_{1}",
    "python {0}/extract_structures_minima.py {1}_{2}_{2}_projection_minima.txt {1}_{2}_common_1.xtc {1}_{2}_common.pdb {1}_{2}_{2}",
]

# Run each command with the prefixes
for cmd in cmds:
    # Format the command with the correct prefixes and the script path
    cmd = cmd.format(script_path, prefix1, prefix2)
    # Use subprocess to execute the command
    subprocess.run(cmd, shell=True, check=True)

