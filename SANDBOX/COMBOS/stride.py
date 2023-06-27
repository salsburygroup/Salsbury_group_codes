import subprocess
import argparse

def run_script(script_path, script_args):
    cmd = ["python", script_path] + script_args
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Error occurred while running the script: {script_path}")
        print(stderr.decode())

    return stdout.decode()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process trajectories.')
    parser.add_argument('prefix', type=str, help='prefix for the files')
    parser.add_argument('n', type=int, help='number of times to run the script')
    parser.add_argument('--pdb_prefix', type=str, default='ionized', help='prefix for the pdb files')
    parser.add_argument('--length', type=int, help='optional length parameter')

    args = parser.parse_args()

    length_10 = 10
    length_100 = 100

    if args.length is not None:
        length_10 = args.length // 100
        length_100 = args.length // 10

    for i in range(1, args.n + 1):
        script_path = "/deac/phy/salsburyGrp/group_python/process_trajectory_juststride.py"
        script_args_10 = [f"{args.prefix}_{i}.xtc", f"{args.pdb_prefix}.pdb", f"{args.prefix}_{length_10}_{i}.xtc", f"{args.prefix}_{length_10}_{i}.pdb", str(length_10), "10000"]
        script_args_100 = [f"{args.prefix}_{i}.xtc", f"{args.pdb_prefix}.pdb", f"{args.prefix}_{length_100}_{i}.xtc", f"{args.prefix}_{length_100}_{i}.pdb", str(length_100), "1000"]

        print(run_script(script_path, script_args_10))
        print(run_script(script_path, script_args_100))

# python stride.py my_prefix 5 --pdb_prefix my_pdb_prefix --length 100
# if pdb_prefix is not passed then ionized.pdb, and if length is not based then 1000
# this code will not run on simulations shorter than 100ns
