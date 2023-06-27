import subprocess
import argparse
import os

def run_script(script_path, script_args):
    cmd = ["python", script_path] + script_args
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Error occurred while running the script: {script_path}")
        print(stderr.decode())

    return stdout.decode()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Wrap trajectories.')
    parser.add_argument('prefix', type=str, help='prefix for the files')
    parser.add_argument('n', type=int, help='number of times to run the script')
    parser.add_argument('--psf_prefix', type=str, default='ionized', help='prefix for the psf files')
    parser.add_argument('--length', type=int, help='optional length parameter')
    parser.add_argument('--path', type=str, default='/deac/phy/salsburyGrp/', help='optional path parameter')

    args = parser.parse_args()

    length_10 = 10
    length_100 = 100

    if args.length is not None:
        length_10 = args.length // 100
        length_100 = args.length // 10

    for i in range(1, args.n + 1):
        script_path = os.path.join(args.path, "group_python/wrap.py")
        script_args_10 = ["-s", f"{args.psf_prefix}.psf", "-t", f"{args.prefix}_{length_10}_{i}.xtc", "-o", f"{args.prefix}_wrapped_{length_10}_{i}.xtc"]
        script_args_100 = ["-s", f"{args.psf_prefix}.psf", "-t", f"{args.prefix}_{length_100}_{i}.xtc", "-o", f"{args.prefix}_wrapped_{length_100}_{i}.xtc"]

        print(run_script(script_path, script_args_10))
        print(run_script(script_path, script_args_100))

