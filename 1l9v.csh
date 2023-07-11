mkdir -p ANALYSIS
mkdir -p ANALYSIS/SCRATCH
mkdir -p ANALYSIS/RMSD
mkdir -p ANALYSIS/TRAJ
mkdir -p ANALYSIS/RMSF
mkdir -p ANALYSIS/CORR
mkdir -p ANALYSIS/PCA
mkdir -p ANALYSIS/START
python /deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py 1l9v_1.xtc ionized.pdb 1l9v_10_1.xtc 1l9v_10_1.pdb 10 1000
/home/salsbufr/anaconda3/bin/conda run -n htmd python /deac/phy/salsburyGrp/group_python/wrap.py -t 1l9v_10_1.xtc -s ionized.psf -o 1l9v_solvated_wrapped_10_1.xtc
python /deac/phy/salsburyGrp/group_python/process_trajectory_large_atomrange.py 1l9v_solvated_wrapped_10_1.xtc 1l9v_10_1.pdb 1l9v_nowat_10_1.xtc 1l9v_nowat_10_1.pdb 1 1000 0-5110
python /deac/phy/salsburyGrp/group_python/calculate_rmsd.py 1l9v_nowat_10_1.xtc 1l9v_nowat_10_1.pdb 1l9v_10_1
python /deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py 1l9v_1.xtc ionized.pdb 1l9v_2_1.xtc 1l9v_2_1.pdb 2 10000
/home/salsbufr/anaconda3/bin/conda run -n htmd python /deac/phy/salsburyGrp/group_python/wrap.py -t 1l9v_2_1.xtc -s ionized.psf -o 1l9v_solvated_wrapped_2_1.xtc
python /deac/phy/salsburyGrp/group_python/process_trajectory_large_atomrange.py 1l9v_solvated_wrapped_2_1.xtc 1l9v_2_1.pdb 1l9v_nowat_2_1.xtc 1l9v_nowat_2_1.pdb 1 10000 0-5110
python /deac/phy/salsburyGrp/group_python/calculate_rmsd.py 1l9v_nowat_2_1.xtc 1l9v_nowat_2_1.pdb 1l9v_2_1
python /deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py 1l9v_2.xtc ionized.pdb 1l9v_10_2.xtc 1l9v_10_2.pdb 10 1000
/home/salsbufr/anaconda3/bin/conda run -n htmd python /deac/phy/salsburyGrp/group_python/wrap.py -t 1l9v_10_2.xtc -s ionized.psf -o 1l9v_solvated_wrapped_10_2.xtc
python /deac/phy/salsburyGrp/group_python/process_trajectory_large_atomrange.py 1l9v_solvated_wrapped_10_2.xtc 1l9v_10_2.pdb 1l9v_nowat_10_2.xtc 1l9v_nowat_10_2.pdb 1 1000 0-5110
python /deac/phy/salsburyGrp/group_python/calculate_rmsd.py 1l9v_nowat_10_2.xtc 1l9v_nowat_10_2.pdb 1l9v_10_2
python /deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py 1l9v_2.xtc ionized.pdb 1l9v_2_2.xtc 1l9v_2_2.pdb 2 10000
/home/salsbufr/anaconda3/bin/conda run -n htmd python /deac/phy/salsburyGrp/group_python/wrap.py -t 1l9v_2_2.xtc -s ionized.psf -o 1l9v_solvated_wrapped_2_2.xtc
python /deac/phy/salsburyGrp/group_python/process_trajectory_large_atomrange.py 1l9v_solvated_wrapped_2_2.xtc 1l9v_2_2.pdb 1l9v_nowat_2_2.xtc 1l9v_nowat_2_2.pdb 1 10000 0-5110
python /deac/phy/salsburyGrp/group_python/calculate_rmsd.py 1l9v_nowat_2_2.xtc 1l9v_nowat_2_2.pdb 1l9v_2_2
python /deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py 1l9v_3.xtc ionized.pdb 1l9v_10_3.xtc 1l9v_10_3.pdb 10 1000
/home/salsbufr/anaconda3/bin/conda run -n htmd python /deac/phy/salsburyGrp/group_python/wrap.py -t 1l9v_10_3.xtc -s ionized.psf -o 1l9v_solvated_wrapped_10_3.xtc
python /deac/phy/salsburyGrp/group_python/process_trajectory_large_atomrange.py 1l9v_solvated_wrapped_10_3.xtc 1l9v_10_3.pdb 1l9v_nowat_10_3.xtc 1l9v_nowat_10_3.pdb 1 1000 0-5110
python /deac/phy/salsburyGrp/group_python/calculate_rmsd.py 1l9v_nowat_10_3.xtc 1l9v_nowat_10_3.pdb 1l9v_10_3
python /deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py 1l9v_3.xtc ionized.pdb 1l9v_2_3.xtc 1l9v_2_3.pdb 2 10000
/home/salsbufr/anaconda3/bin/conda run -n htmd python /deac/phy/salsburyGrp/group_python/wrap.py -t 1l9v_2_3.xtc -s ionized.psf -o 1l9v_solvated_wrapped_2_3.xtc
python /deac/phy/salsburyGrp/group_python/process_trajectory_large_atomrange.py 1l9v_solvated_wrapped_2_3.xtc 1l9v_2_3.pdb 1l9v_nowat_2_3.xtc 1l9v_nowat_2_3.pdb 1 10000 0-5110
python /deac/phy/salsburyGrp/group_python/calculate_rmsd.py 1l9v_nowat_2_3.xtc 1l9v_nowat_2_3.pdb 1l9v_2_3
python /deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py 1l9v_4.xtc ionized.pdb 1l9v_10_4.xtc 1l9v_10_4.pdb 10 1000
/home/salsbufr/anaconda3/bin/conda run -n htmd python /deac/phy/salsburyGrp/group_python/wrap.py -t 1l9v_10_4.xtc -s ionized.psf -o 1l9v_solvated_wrapped_10_4.xtc
python /deac/phy/salsburyGrp/group_python/process_trajectory_large_atomrange.py 1l9v_solvated_wrapped_10_4.xtc 1l9v_10_4.pdb 1l9v_nowat_10_4.xtc 1l9v_nowat_10_4.pdb 1 1000 0-5110
python /deac/phy/salsburyGrp/group_python/calculate_rmsd.py 1l9v_nowat_10_4.xtc 1l9v_nowat_10_4.pdb 1l9v_10_4
python /deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py 1l9v_4.xtc ionized.pdb 1l9v_2_4.xtc 1l9v_2_4.pdb 2 10000
/home/salsbufr/anaconda3/bin/conda run -n htmd python /deac/phy/salsburyGrp/group_python/wrap.py -t 1l9v_2_4.xtc -s ionized.psf -o 1l9v_solvated_wrapped_2_4.xtc
python /deac/phy/salsburyGrp/group_python/process_trajectory_large_atomrange.py 1l9v_solvated_wrapped_2_4.xtc 1l9v_2_4.pdb 1l9v_nowat_2_4.xtc 1l9v_nowat_2_4.pdb 1 10000 0-5110
python /deac/phy/salsburyGrp/group_python/calculate_rmsd.py 1l9v_nowat_2_4.xtc 1l9v_nowat_2_4.pdb 1l9v_2_4
python /deac/phy/salsburyGrp/group_python/merge_trajectories.py 1l9v_nowat_10_1.xtc 1l9v_nowat_10_2.xtc 1l9v_nowat_10_3.xtc 1l9v_nowat_10_4.xtc --topology 1l9v_autopsf.psf --reference_structure 1l9v_nowat_10_3.pdb --output 1l9v_nowat_10 --chunk_size 1000 --selection "all"
python /deac/phy/salsburyGrp/group_python/merge_trajectories.py 1l9v_nowat_2_1.xtc 1l9v_nowat_2_2.xtc 1l9v_nowat_2_3.xtc 1l9v_nowat_2_4.xtc --topology 1l9v_autopsf.psf --reference_structure 1l9v_nowat_2_3.pdb --output 1l9v_nowat_2 --chunk_size 10000 --selection "all"
mv *_rmsd.txt ANALYSIS/RMSD/
mv *_rmsd.png ANALYSIS/RMSD/
mv 1l9v_nowat_10_merged.xtc ANALYSIS/TRAJ/
mv 1l9v_nowat_2_merged.xtc ANALYSIS/TRAJ/
mv 1l9v_nowat_2_1.pdb ANALYSIS/TRAJ/
mv 1l9v_1.xtc ANALYSIS/START/
mv 1l9v_2.xtc ANALYSIS/START/
mv 1l9v_3.xtc ANALYSIS/START/
mv 1l9v_4.xtc ANALYSIS/START/
mv ionized.pdb ANALYSIS/START/
mv ionized.psf ANALYSIS/START/
mv 1l9v_autopsf.psf ANALYSIS/START/
mv *.xtc ANALYSIS/SCRATCH/
mv *.pdb ANALYSIS/SCRATCH/
python /deac/phy/salsburyGrp/group_python/calculate_rmsf.py ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb ANALYSIS/TRAJ/1l9v_nowat_2_merged.xtc 1l9v
mv *rmsf* ANALYSIS/RMSF/
python /deac/phy/salsburyGrp/group_python/calculate_corr.py  ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb ANALYSIS/TRAJ/1l9v_nowat_2_merged.xtc  1l9v_2 --align 
mv *correlation_matrix* ANALYSIS/CORR/
python /deac/phy/salsburyGrp/group_python/calculate_pca_projections.py ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb ANALYSIS/TRAJ/1l9v_nowat_2_merged.xtc all 1l9v
python /deac/phy/salsburyGrp/group_python/calculate_pca_corner.py ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb ANALYSIS/TRAJ/1l9v_nowat_2_merged.xtc all 1l9v 5
python /deac/phy/salsburyGrp/group_python/calculate_FES_PCA.py  1l9v_projections.txt 1l9v -b 15
python /deac/phy/salsburyGrp/group_python/find_minima_find_structures.py 1l9v_bins_15_free_energy.txt 1l9v_bin_indices.txt 1l9v
python /deac/phy/salsburyGrp/group_python/extract_structures_minima.py 1l9v_projection_minima.txt ANALYSIS/TRAJ/1l9v_nowat_2_merged.xtc ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb 1l9v
mv *projection* ANALYSIS/PCA/
mv *energy* ANALYSIS/PCA/
mv *minima* ANALYSIS/PCA/
mv *bin* ANALYSIS/PCA/
mv *corner* ANALYSIS/PCA/
python /deac/phy/salsburyGrp/group_python/cluster_HDBSCAN.py ANALYSIS/TRAJ/1l9v_nowat_10_merged.xtc ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb 1l9v
python /deac/phy/salsburyGrp/group_python/cluster_traj.py ANALYSIS/TRAJ/1l9v_nowat_10_merged.xtc ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb 1l9v 10 1l9v_HDBSCAN_frame_clusters.txt
mkdir -p ANALYSIS/HDBSCAN
mv *HDBSCAN* ANALYSIS/HDBSCAN/
mv 1l9v/* ANALYSIS/HDBSCAN/
rmdir 1l9v
python /deac/phy/salsburyGrp/group_python/cluster_AH.py ANALYSIS/TRAJ/1l9v_nowat_10_merged.xtc ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb 1l9v
python /deac/phy/salsburyGrp/group_python/cluster_traj.py ANALYSIS/TRAJ/1l9v_nowat_10_merged.xtc ANALYSIS/TRAJ/1l9v_nowat_2_1.pdb 1l9v 10 1l9v_AH_frame_clusters.txt
mkdir -p ANALYSIS/AH
mv *AH* ANALYSIS/AH/
mv 1l9v/* ANALYSIS/AH/
rmdir 1l9v
python /deac/phy/salsburyGrp/group_python/HBONDS/find_hbonds.py ANALYSIS/START/1l9v_autopsf.psf  ANALYSIS/TRAJ/1l9v_nowat_2_merged.xtc 1l9v
python /deac/phy/salsburyGrp/group_python/HBONDS/filter_hbonds_cutoff.py 1l9v
python /deac/phy/salsburyGrp/group_python/HBONDS/cluster_hbonds.py 1l9v
mkdir -p ANALYSIS/HBONDS
mv *hbond* ANALYSIS/HBONDS/
mv *info* ANALYSIS/HBONDS/
