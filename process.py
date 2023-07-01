import os
import argparse
import glob
import math
import subprocess

def main(prefix, n, atom_range, path, pdb_prefix, psf_prefix, length, nowat_psf_prefix, merge_selection, norun, corr_range, pca_selection, bins, conda_path, n_clusters):

    first_value = length // 10
    second_value = length // 100
    first_multiplier = 1000
    second_multiplier = 10000
    smaller_number = min(first_value, second_value)
    larger_number = max(first_value, second_value)
    script_path = os.path.join(path, "group_python")
    analysis_dir = f"ANALYSIS"
    scratch_dir = os.path.join(analysis_dir, "SCRATCH")
    rmsd_dir = os.path.join(analysis_dir, "RMSD")
    traj_dir = os.path.join(analysis_dir, "TRAJ")
    rmsf_dir = os.path.join(analysis_dir, "RMSF")
    corr_dir = os.path.join(analysis_dir, "CORR")
    pca_dir = os.path.join(analysis_dir, "PCA")
    start_dir = os.path.join(analysis_dir, "START")
    commands = []

    # Directory creation commands
    dir_commands = [f"mkdir -p {directory}" for directory in [analysis_dir, scratch_dir, rmsd_dir, traj_dir, rmsf_dir, corr_dir, pca_dir, start_dir]]
    commands.extend(dir_commands)

    if not args.analysisonly:
     for i in range(1, n + 1):
        for value, multiplier in zip([first_value, second_value], [first_multiplier, second_multiplier]):
            # process_trajectory_juststride.py
            command = f"python {os.path.join(script_path, 'process_trajectory_juststride.py')} {prefix}_{i}.xtc {pdb_prefix}.pdb {prefix}_{value}_{i}.xtc {prefix}_{value}_{i}.pdb {value} {multiplier}"
            commands.append(command)

            # wrap.py
            htmd_conda_path = "/home/salsbufr/anaconda3/bin/"  # Replace with your path to conda's bin
            command = f"{conda_path}conda run -n htmd python {os.path.join(script_path, 'wrap.py')} -t {prefix}_{value}_{i}.xtc -s {psf_prefix}.psf -o {prefix}_solvated_wrapped_{value}_{i}.xtc"
            commands.append(command)

            # process_trajectory_large_atomrange.py
            command = f"python {os.path.join(script_path, 'process_trajectory_large_atomrange.py')} {prefix}_solvated_wrapped_{value}_{i}.xtc {prefix}_{value}_{i}.pdb {prefix}_nowat_{value}_{i}.xtc {prefix}_nowat_{value}_{i}.pdb 1 {multiplier} {atom_range}"
            commands.append(command)

            # calculate_rmsd.py
            command = f"python {os.path.join(script_path, 'calculate_rmsd.py')} {prefix}_nowat_{value}_{i}.xtc {prefix}_nowat_{value}_{i}.pdb {prefix}_{value}_{i}"
            commands.append(command)

    # merge_trajectories.py
     for value, multiplier in zip([first_value, second_value], [first_multiplier, second_multiplier]):
        input_xtcs = " ".join(f"{prefix}_nowat_{value}_{i}.xtc" for i in range(1, n + 1))
        command = f"python {os.path.join(script_path, 'merge_trajectories.py')} {input_xtcs} --topology {nowat_psf_prefix}.psf --reference_structure {prefix}_nowat_{value}_3.pdb --output {prefix}_nowat_{value} --chunk_size {multiplier} --selection \"{args.merge_selection}\""
        commands.append(command)

     # File moving commands
     # Move RMSD files to RMSD directory
     commands.append(f"mv *_rmsd.txt {rmsd_dir}/")
     commands.append(f"mv *_rmsd.png {rmsd_dir}/")

     # Move two merged xtc files and a reference pdb file to TRAJ directory
     commands.append(f"mv {prefix}_nowat_{first_value}_merged.xtc {traj_dir}/")
     commands.append(f"mv {prefix}_nowat_{second_value}_merged.xtc {traj_dir}/")
     commands.append(f"mv {prefix}_nowat_{smaller_number}_1.pdb {traj_dir}/")
     for i in range(1, n + 1):
       commands.append(f"mv {prefix}_{i}.xtc {start_dir}/")

     # Move pdb and psf files to START directory
     commands.append(f"mv {pdb_prefix}.pdb {start_dir}/")
     commands.append(f"mv {psf_prefix}.psf {start_dir}/")
     commands.append(f"mv {nowat_psf_prefix}.psf {start_dir}/")
 
     # Move remaining files to SCRATCH directory
     commands.append(f"mv *.xtc {scratch_dir}/")
     commands.append(f"mv *.pdb {scratch_dir}/")

     # calculate_rmsf.py
    command = f"python {os.path.join(script_path, 'calculate_rmsf.py')} {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc {prefix}"
    commands.append(command)
   
    #Move RMSF files to RMSF directory
    commands.append(f"mv *rmsf* {rmsf_dir}/")

    # calculate_corr.py
    if corr_range is not None:
      command = f"python {os.path.join(script_path, 'calculate_corr.py')} {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb  {traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc {prefix}_{smaller_number} --range {corr_range}"
    else:
      command = f"python {os.path.join(script_path, 'calculate_corr.py')}  {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc  {prefix}_{smaller_number} --align "
    commands.append(command)

    # Move correlation matrix files to CORR directory
    commands.append(f"mv *correlation_matrix* {corr_dir}/")

    # calculate_pca_projections.py
    command = f"python {os.path.join(script_path, 'calculate_pca_projections.py')} {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc {pca_selection} {prefix}"
    commands.append(command)

   # corner plot
    command = f"python {os.path.join(script_path, 'calculate_pca_corner.py')} {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc {pca_selection} {prefix} 5"
    commands.append(command)
 
   # Calculate bin number
    bin_number = round(1 + math.log2(n*length*100/smaller_number))
    bins = bins if bins is not None else bin_number

   # calculate_FES_PCA.py
    command = f"python {os.path.join(script_path, 'calculate_FES_PCA.py')}  {prefix}_projections.txt {prefix} -b {bins}"
    commands.append(command)

   # find_minima_find_structures.py
    command = f"python {os.path.join(script_path, 'find_minima_find_structures.py')} {prefix}_bins_{bins}_free_energy.txt {prefix}_bin_indices.txt {prefix}"
    commands.append(command)

    # extract_structures_minima.py
    command = f"python {os.path.join(script_path, 'extract_structures_minima.py')} {prefix}_projection_minima.txt {traj_dir}/{prefix}_nowat_{smaller_number}_merged.xtc {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {prefix}"
    commands.append(command)

    # Move PCA files to PCA directory
    commands.append(f"mv *projection* {pca_dir}/")
    commands.append(f"mv *energy* {pca_dir}/")
    commands.append(f"mv *minima* {pca_dir}/")
    commands.append(f"mv *bin* {pca_dir}/")
    commands.append(f"mv *corner* {pca_dir}/") 
    commands.append(f"mv *contour* {pca_dir}/")
    # HDBSCAN Clustering
    command = f"python {os.path.join(script_path, 'cluster_HDBSCAN.py')} {traj_dir}/{prefix}_nowat_{larger_number}_merged.xtc {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {prefix}"
    commands.append(command)

    # HDBSCAN Cluster Trajectory Analysis
    command = f"python {os.path.join(script_path, 'cluster_traj.py')} {traj_dir}/{prefix}_nowat_{larger_number}_merged.xtc {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {prefix} {n_clusters} {prefix}_HDBSCAN_frame_clusters.txt"
    commands.append(command)

    # HDBSCAN directory creation command
    hdbscan_dir = os.path.join(analysis_dir, "HDBSCAN")
    commands.append(f"mkdir -p {hdbscan_dir}")

   # Move files with HDBSCAN in their names into the HDBSCAN directory
    commands.append(f"mv *HDBSCAN* {hdbscan_dir}/")

    # Move everything from the prefix directory into the HDBSCAN directory and remove the prefix directory
    commands.append(f"mv {prefix}/* {hdbscan_dir}/")
    commands.append(f"rmdir {prefix}")

    # AH Clustering
    command = f"python {os.path.join(script_path, 'cluster_AH.py')} {traj_dir}/{prefix}_nowat_{larger_number}_merged.xtc {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {prefix}"
    commands.append(command)

    # AH Cluster Trajectory Analysis
    command = f"python {os.path.join(script_path, 'cluster_traj.py')} {traj_dir}/{prefix}_nowat_{larger_number}_merged.xtc {traj_dir}/{prefix}_nowat_{smaller_number}_1.pdb {prefix} {n_clusters} {prefix}_AH_frame_clusters.txt"
    commands.append(command)

    # AH directory creation command
    AH_dir = os.path.join(analysis_dir, "AH")
    commands.append(f"mkdir -p {AH_dir}")

    # Move files with AH in their names into the AH directory
    commands.append(f"mv *AH* {AH_dir}/")

    # Move everything from the prefix directory into the AH directory and remove the prefix directory
    commands.append(f"mv {prefix}/* {AH_dir}/")
    commands.append(f"rmdir {prefix}")
 
    if norun:
        with open(f"{prefix}.csh", 'w') as f:
            for command in commands:
                f.write(command + "\n")
    else:
        for command in dir_commands:
            os.system(command)

        for command in commands:
            os.system(command)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to perform MD analysis.")
    parser.add_argument('prefix', type=str, help='Prefix for the simulation.')
    parser.add_argument('-n', '--n', type=int, default=4, help='Number of chunks for the simulation. Default is 4.')
    parser.add_argument('atom_range', type=str, help='Range of atoms to consider in the simulation.')
    parser.add_argument('--path', type=str, default="/deac/phy/salsburyGrp/", help='Path to group_python scripts.')
    parser.add_argument('--pdb_prefix', type=str, default="ionized", help='Prefix for the pdb files.')
    parser.add_argument('--psf_prefix', type=str, default="ionized", help='Prefix for the psf files.')
    parser.add_argument('--length', type=int, default=1000, help='Total length of the trajectory.')
    # Parse the known arguments to get the prefix value
    known_args, remaining = parser.parse_known_args()
    # Add the '--nowat_psf_prefix' argument using the prefix value
    parser.add_argument('--nowat_psf_prefix', type=str, default=f"{known_args.prefix}_autopsf", help='Prefix for the no-water psf files.')
    parser.add_argument('--merge_selection', type=str, default="all", help='Selection for merge_trajectories.py script.')
    parser.add_argument('--norun', action='store_true', help='If provided, the script only writes the commands to a .csh file instead of running them.')
    parser.add_argument('--corr_range', type=str, help='Range for calculate_corr.py script.')
    parser.add_argument('--pca_selection', type=str, default="all", help='Selection for calculate_pca_projections.py script.')
    parser.add_argument('--bins', type=int, help='Number of bins for calculate_FES_PCA.py and find_minima_find_structures.py scripts.')
    parser.add_argument('--conda_path', type=str, default=os.path.join('/home', os.getlogin(), 'anaconda3/bin/'), help='Path to conda\'s bin.')
    parser.add_argument('--analysisonly', action='store_true', help='If provided, the script only runs the analysis part.')
    parser.add_argument("--n_clusters", type=int, default=10)
    args = parser.parse_args()

if __name__ == "__main__":
    args = parser.parse_args()

main(args.prefix, args.n, args.atom_range, args.path, args.pdb_prefix, args.psf_prefix, args.length, args.nowat_psf_prefix, args.merge_selection, args.norun, args.corr_range, args.pca_selection, args.bins, args.conda_path, args.n_clusters)

