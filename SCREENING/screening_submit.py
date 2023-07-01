import os
import time
import math
import sys

def submit_jobs(screen_type):
    # determine the number of lines in the file
    with open(f'/deac/phy/salsburyGrp/autodock/ZINC_libraries_06232023/ZINC_clean_instock_{screen_type}like.smi', 'r') as f:
        lines = f.readlines()
        library_size = len(lines)

    # calculate the ceiling integer of library_size / 1000
    num_jobs = math.ceil(library_size / 100000)

    # loop over the range and submit the jobs
    for n in range(1, num_jobs+1):
        os.system(f'sbatch --export=n={n},m=1000 docking_{screen_type}.slurm')
        time.sleep(2)  # wait for 2 seconds

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python screening_submit.py <screen_type>")
        sys.exit(1)

    screen_type = sys.argv[1]
    submit_jobs(screen_type)

