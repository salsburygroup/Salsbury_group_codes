import numpy as np
import argparse
import subprocess
import os

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description = 'Runs cluster trials over variety of methods and metrics', add_help=False) 

#List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel', action='store', dest='sel',help='atoms',type=str,required=True)
inputs.add_argument('-title', action='store', dest='title',help='Title',type=str,required=True)
inputs.add_argument('-tm', action='store', dest='time', help='Time interval(ps)', type=str, default='10')
inputs.add_argument('-nm', action='store', dest='nm',help='Output file name',type=str,required=True)
inputs.add_argument('-o', action='store', dest='out_dir',help='Output directory',type=str,required=True)

#Parse into useful form
UserInput=parser.parse_args()

# Find Helper scripts
cluster_script = '/home/wud18/python/SodiumLoopIonOnOff.py'

# Calculating distance
if not os.path.exists(UserInput.out_dir):
    os.makedirs(UserInput.out_dir)
python_NaBinding_cmd = (
    'sbatch --export=t=' + UserInput.trajectory +
    ',str=' + UserInput.structure +
    ',s=' + cluster_script +
    ',a=' + UserInput.sel + 
    ',b=' + UserInput.title + 
    ',tm=' + UserInput.time +
    ',nm=' + UserInput.nm +
    ',o=' + UserInput.out_dir + 
    ' /home/wud18/bash/NaBindingSubmit.slurm'
    )
print(python_NaBinding_cmd)
subprocess.call(python_NaBinding_cmd, shell=True)
