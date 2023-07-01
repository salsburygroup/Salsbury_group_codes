import sys
import os

def generate_scripts(screen_type, center, pdb_prefix, project_name):
    centers = center.split()

    # Generate ligand_prep.sh
    with open('ligand_prep.sh', 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('for((i=$1;i<=$1+$2;i++))\n')
        f.write('do\n')
        f.write(f'  cd /deac/phy/salsburyGrp/autodock/{project_name}/screen_{screen_type}/screen_$n/zinc\"$i\"\n')
        f.write('/home/luy/MGLTools-1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l zinc.pdb\n')
        f.write('done\n')

    # Generate run_vina.sh
    with open('run_vina.sh', 'w') as f:
        f.write('#!/bin/bash -l\n')
        f.write('. /etc/profile.d/modules.sh\n')
        f.write('module load rhel7/vina/1.1.2-intel-2018\n')
        f.write('for((i=$1;i<=$1+$2;i++))\n')
        f.write('do\n')
        f.write(f'cd /deac/phy/salsburyGrp/autodock/{project_name}/screen_{screen_type}/screen_$n/zinc\"$i\"\n')
        f.write(f'vina --receptor {pdb_prefix}.pdbqt --ligand zinc.pdbqt --center_x {centers[0]} --center_y {centers[1]}  --center_z {centers[2]} --size_x 30 --size_y 30 --size_z 30 --cpu 1 > log_zinc\"$i\"\n')
        f.write('done\n')

    # Generate docking_{screen}.slurm
    with open(f'docking_{screen_type}.slurm', 'w') as f:
        f.write('#!/bin/bash -l\n')
        f.write('#SBATCH --partition=medium\n')
        f.write('#SBATCH --job-name="docking_DT.job"\n')
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --tasks-per-node=1\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write('#SBATCH --time=7-00:00:00\n')
        f.write('#SBATCH --account="salsburyGrp"\n')
        f.write('##SBATCH --mail-type=BEGIN,END,FAIL\n')
        f.write('##SBATCH --mail-user=salsbufr@wfu.edu\n')
        f.write('#SBATCH --output="SLURM_OUTS/docking_DT_lead_core_1_%j.o"\n')
        f.write('#SBATCH --error="SLURM_ERRORS/docking_DT_lead_core_1_%j.e"\n')
        f.write('#SBATCH --mem=1gb\n')
        f.write('source /etc/profile\n')
        f.write(f'mkdir /deac/phy/salsburyGrp/autodock/{project_name}\n')
        f.write(f'mkdir /deac/phy/salsburyGrp/autodock/{project_name}/screen_{screen_type}\n')
        f.write(f'cd /deac/phy/salsburyGrp/autodock/{project_name}/screen_{screen_type}\n')
        f.write('mkdir RESULTS\n')
        f.write('mkdir screen_$n\n')
        f.write('cd screen_$n\n')
        f.write('for((i=$n;i<=$n+$m;i++))\n')
        f.write('do\n')
        f.write('mkdir zinc"$i"\n')
        f.write('cd zinc"$i"\n')
        f.write('rm -r *\n')
        f.write(f'cp /home/salsbufr/{project_name}/{pdb_prefix}.pdbqt .\n')
        f.write('/home/salsbufr/babel/bin/babel -i smi /deac/phy/salsburyGrp/autodock/ZINC_libraries_06232023/ZINC_clean_instock_' + screen_type + 'like.smi -O zinc.pdb  -f $i -l $i --gen3d\n')

        f.write(f'sh /home/salsbufr/{project_name}/SCREENING/ligand_prep.sh $i 0\n')
        f.write(f'sh /home/salsbufr/{project_name}/SCREENING/run_vina.sh $i 0\n')
        f.write(f'awk \'{{if (NF==4 && $2<-8.2 ){{print $2,$1,FILENAME}}}}\' log_* >>/deac/phy/salsburyGrp/autodock/{project_name}/screen_{screen_type}/best_hits.dat\n')
        f.write('awk \'{if (NF==4 && $2<-8.2 ){print $2,$1,FILENAME}}\' log_* >> hits.dat\n')
        f.write('if [ -s hits.dat ]\n')
        f.write('then\n')
        f.write('sed \'/ENDROOT/d\' zinc_out.pdbqt > tmp\n')
        f.write('sed \'/ENDBRANCH/d\' tmp > tmp2\n')
        f.write(f'cp tmp2 /deac/phy/salsburyGrp/autodock/{project_name}/screen_{screen_type}/RESULTS/zinc_out_${{i}}.pdb\n')
        f.write('fi\n')
        f.write(f'cd /deac/phy/salsburyGrp/autodock/{project_name}/screen_{screen_type}/screen_$n\n')
        f.write('rm -r zinc"$i"\n')
        f.write('done\n')

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <screen_type> <center> <pdb_prefix> <project_name>")
        sys.exit(1)

    screen_type = sys.argv[1]
    center = sys.argv[2]
    pdb_prefix = sys.argv[3]
    project_name = sys.argv[4]

generate_scripts(screen_type, center, pdb_prefix, project_name)
