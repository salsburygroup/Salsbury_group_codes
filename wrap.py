from htmd.ui import *
import argparse
#htmd_register()
parser = argparse.ArgumentParser(description='Calculate, save and plot RMSF', add_help=False)
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s',
                    action='store',
                    dest='structure',
                    help='Structure file corresponding to trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-t',
                    action='store',
                    dest='trajectory',
                    help='Trajectory',
                    type=str,
                    required=True)

inputs.add_argument('-sel',
                    action='store',
                    dest='sel',
                    help='Atom selection',
                    type=str,
                    default='protein')

inputs.add_argument('-o',
                    action='store',
                    dest='out_dir',
                    help='Output prefix for text and png',
                    type=str,
                    required=True)
UserInput = parser.parse_args()
mol = Molecule(UserInput.structure)
mol.read(UserInput.trajectory)
mol.wrap(UserInput.sel) # or "protein and resid X to Y" if you want to wrap around specific protein residues
mol.write(UserInput.out_dir)



