#! /usr/bin/env python

#####
#hbonds.py
#Thu Apr 23 15:43:40 EDT 2015
#Ryan Melvin
#####
#Credit:https://mdanalysis.googlecode.com/svn/trunk/doc/html/documentation_pages/analysis/hbonds.html
#If you use this script, cite
#N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319-2327. doi:10.1002/jcc.21787
#AND
#L.M. Gregoret, S.D. Rader, R.J. Fletterick, and F.E. Cohen. Hydrogen bonds involving sulfur atoms in proteins. Proteins, 9(2):99-107, 1991. 10.1002/prot.340090204.
####
#TODO:
#####
#Example call
#python /Users/melvrl13/Documents/AMD/AMD-PYTHON/Analysis/hbond_analysis.py -structure /Volumes/RyanMdata/sufCandD/RelaxationSimulations/SufCD_with_ATP/Homology/LastFrameOfHomologyRound2/Complex.psf -t /Volumes/RyanMdata/sufCandD/RelaxationSimulations/SufCD_with_ATP/Homology/LastFrameOfHomologyRound2/round1/SufCD_ATP_relax_strip_stride.dcd --matchVMD -sel1 all -sel2 all -sel1_type both -d 4 -a 60 -o /Volumes/RyanMdata/sufCandD/RelaxationSimulations/SufCD_with_ATP/Homology/LastFrameOfHomologyRound2/round1/hbond_test.txt
#MDANALYSIS default atom list https://mdanalysis.googlecode.com/svn/trunk/doc/html/documentation_pages/analysis/hbonds.html

#Outputs: Hydrogen bonds table with columns 0)time 1)donor index 2) acceptor index 3) donor residue name 4) donor residue id 5) donor atom 6) acceptor residue name 7) acceptor residue id 8) acceptor atom 9) distance 10) angle
#NOTE 1-based atom indexing
#NOTE 0-based column indexing

# Dependencies
from __future__ import division
import MDAnalysis
import MDAnalysis.analysis.hbonds
import argparse
import sys
import pandas as pd
import numpy as np
from warnings import warn

print("Please cite: "\
    "N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: "\
    "A Toolkit for the Analysis of Molecular Dynamics Simulations. "\
    "J. Comput. Chem. 32 (2011), 2319-2327. doi:10.1002/jcc.21787 "\
    "AND "\
    "L.M. Gregoret, S.D. Rader, R.J. Fletterick, and F.E. Cohen. "\
    "Hydrogen bonds involving sulfur atoms in proteins. Proteins, 9(2):99-107, 1991. 10.1002/prot.340090204.")

#Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(
        description = (
            'output columns: 0)time 1)donor index 2) acceptor index 3) donor residue name 4) donor residue id 5) donor atom 6) acceptor residue name 7) acceptor residue id 8) acceptor atom 9) distance 10) angle'
            ), 
        add_help=False
        ) 

# List all possible user input
inputs=parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-s', action='store', dest='structure',help='Structure file corresponding to trajectory',type=str,required=True)
inputs.add_argument('-t', action='store', dest='trajectory',help='Trajectory',type=str,required=True)
inputs.add_argument('-sel1', action='store', dest='sel1', help='Atom selection 1',type=str,default='protein')
inputs.add_argument('-sel2', action='store', dest='sel2', help='Atom selection 2',type=str,default='all')
inputs.add_argument('-sel1_type', action='store', dest='sel1_type', help='Selection 1 type, i.e. donor, acceptor or both',type=str,default='both')
inputs.add_argument('-update_sel1', action='store_true', dest='update_sel1', help='Upate selection 1 each frame?',default=False)
inputs.add_argument('-update_sel2', action='store_true', dest='update_sel2', help='Upate selection 2 each frame?',default=False)
inputs.add_argument('-d', '--distance', action='store', dest='distance',help='Distance cutoff in angstroms',type=float,default=3.2)
inputs.add_argument('-a', '--angle',  action='store', dest='angle', help='Angle cutoff in degrees', type=float, default=120.00)
inputs.add_argument('-o', action='store', dest='out_name',help='Output prefix with no extension',type=str,required=True)
inputs.add_argument('--extra_acceptors', action='store', dest='acceptors',help='Acceptors in addition to those defined by charmm27',type=str,default=None)
inputs.add_argument('--extra_donors', action='store', dest='donors',help='Donors in addition to those defined by charmm27',type=str,default=None)
inputs.add_argument('--distance_type', action='store', dest='distance_type',help='heavy or hydrogen?',type=str,default='heavy')
inputs.add_argument('--polar', action='store_true', help='Emulates VMD polar atom only algorithm. Ignores extra donors and acceptors options. Also, This options will not work for atom names longer than three characters.')
# Three letter atom name restriction to be fixed in future version.



# Parse into useful form
UserInput=parser.parse_args()

# Define the universe (i.e., molecule in VMD)
u = MDAnalysis.Universe(UserInput.structure, UserInput.trajectory, permissive=True)

# Emulate VMD's detection algorithm if option is selected
if UserInput.polar:
    warn("The polar flag ignores additional donors or acceptors. Instead, it searches all polar atoms -- N, O, S and F.")

    donors_acceptors = u.select_atoms("name O or name O* or name O** or name N or name N* or name N** or name S or name S* or name S** or name F or name F* or name F**")
    donors_acceptors_list = np.unique(donors_acceptors.names)

    # Setup MD Analysis input
    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(
        u, 
        UserInput.sel1, 
        UserInput.sel2, 
        selection1_type=UserInput.sel1_type, 
        update_selection1=UserInput.update_sel1, 
        update_selection2=UserInput.update_sel2, 
        distance=UserInput.distance, 
        angle=UserInput.angle, 
        donors=donors_acceptors_list,
        acceptors=donors_acceptors_list,
        detect_hydrogens='distance',
        distance_type=UserInput.distance_type
        )

else:
    # Setup MD Analysis input
    h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(
        u, 
        UserInput.sel1, 
        UserInput.sel2, 
        selection1_type=UserInput.sel1_type, 
        update_selection1=UserInput.update_sel1, 
        update_selection2=UserInput.update_sel2, 
        distance=UserInput.distance, 
        angle=UserInput.angle, 
        donors=UserInput.donors,
        acceptors=UserInput.acceptors,
        detect_hydrogens='distance',
        distance_type=UserInput.distance_type
        )



# Execute analysis.
h.run()

# Format as a table
h.generate_table()

# Save as csv and change tiem to frames
df = pd.DataFrame.from_records(h.table)
df['time']=df['time'].apply(lambda x: int(round(x/u.trajectory.dt))) #Convert from MDAnalysis default time values to frames
df.to_csv(UserInput.out_name + '_raw.csv',index=False) #Save as csv


# Repeat for each frame individually and record
hbond_trajectory = pd.DataFrame(index=list(range(0, len(u.trajectory)))) #Empty data frame in memory
for frame in list(range(0, len(u.trajectory))):
    current_frame_hbonds = df.loc[df['time'] == int(frame), ['donor_resnm', 'donor_resid', 'acceptor_resnm', 'acceptor_resid']] #All hbonds in frame
    current_frame_hbond_pairs = [row['donor_resnm'] + str(row['donor_resid']) + '-' + row['acceptor_resnm'] + str(row['acceptor_resid'])
                                 for index, row in  current_frame_hbonds.iterrows()]
    for pair in current_frame_hbond_pairs:
        hbond_trajectory.loc[frame, pair] = 1 #Fill in the empties with 1 if the hbond occurs
    progress = "\r Motif calculation on Frame " + str(frame) + " of " + str(len(u.trajectory)) #status
    sys.stdout.write(progress)
    sys.stdout.flush() #report status to terminal output
hbond_trajectory = hbond_trajectory.fillna(0) #Fill any leftover empties with 0
hbond_trajectory.to_csv(UserInput.out_name + '_trajectory.csv',index=False) #Save as csv
