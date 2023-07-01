import os
import subprocess
import sys

# Check if directory argument was passed
if len(sys.argv) < 2:
    print("Usage: python script.py /path/to/directory")
    sys.exit(1)

# Initialize total CPU time
total_cpu_time = 0

# Directory containing the output files
directory = sys.argv[1]

# Costs per cpu-hour
direct_cost_rate = 0.03
indirect_cost_rate = 0.02
total_cost_rate = 0.05

# Loop over all .out files in the directory starting with 'slurm'
for filename in os.listdir(directory):
    if filename.startswith("slurm") and filename.endswith(".out"):
        # Extract job ID from the filename
        job_id = filename[5:-4]
        if job_id.startswith("-"):
            job_id = job_id[1:]

        # Get the total CPU time for this job
        sacct_command = f"sacct -j {job_id} --format=TotalCPU --noheader"
        sacct_output = subprocess.check_output(sacct_command, shell=True).decode().strip().split('\n')[0]
        hours, minutes, seconds = map(int, sacct_output.split(':'))

        # Convert CPU time to seconds and add to total
        cpu_time_in_seconds = hours * 3600 + minutes * 60 + seconds
        total_cpu_time += cpu_time_in_seconds

# Convert total CPU time to hours and round it to the nearest second decimal point
total_cpu_time_in_hours = round(total_cpu_time / 3600, 2)

# Calculate costs and round them to the nearest second decimal point
total_cost = round(total_cpu_time_in_hours * total_cost_rate, 2)
direct_cost = round(total_cpu_time_in_hours * direct_cost_rate, 2)
indirect_cost = round(total_cpu_time_in_hours * indirect_cost_rate, 2)

# Print total CPU time and costs
print(f"Total CPU time: {total_cpu_time} seconds, or {total_cpu_time_in_hours} hours")
print(f"Total cost: ${total_cost}")
print(f"Direct cost: ${direct_cost}")
print(f"Indirect cost: ${indirect_cost}")

