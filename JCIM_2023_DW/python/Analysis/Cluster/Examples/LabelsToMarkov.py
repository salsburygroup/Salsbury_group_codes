#!/usr/env python
"""
An example script that takes labels from previous clustering and generates an MSM using pyEmma.
Example command
python LabelsToMarkov.py -l RLM/QT_labels.txt -o /Users/melvrl13/Documents/groupStuff/Jiajie/MelvinQT/RLM
"""


from Analysis.Cluster import Markov, Plotter
import argparse
import numpy
import os
import pyemma.plots as mplt
import matplotlib.pyplot as plt


#def main():
    # Initialize parser. The default help has poor labeling. See http://bugs.python.org/issue9694
parser = argparse.ArgumentParser(description='Convert labels to Markov Rate Matrix', add_help=False)
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-l', action='store', dest='labels_file', help='File containing labels ordered by frame',
                        type=str, required=True)
inputs.add_argument('-o', action='store', dest='output', help='Directory to save output', required=True)
user_input = parser.parse_args()

labels = numpy.genfromtxt(user_input.labels_file, dtype='int')
if -1 in labels:
    labels[labels < 0] = labels.max() + 1

msm = Markov.Model(labels=labels).estimate()

rate_matrix_image_file = os.path.join(user_input.output, 'rate_matrix.png')
Plotter.RateMatrix(out_name=rate_matrix_image_file, msm=msm).plot()

tpt_file = os.path.join(user_input.output, 'transition_path.png')
Plotter.TransitionPath(out_name=tpt_file, msm=msm).plot()

rate_matrix_text_file = os.path.join(user_input.output, 'rate_matrix.txt')
numpy.savetxt(fname=rate_matrix_text_file, X=msm.transition_matrix)

cmtest_file = os.path.join(user_input.output, 'chapman_kolmogorov.png')
cktest = msm.cktest(2)
mplt.plot_cktest(cktest=cktest)
plt.savefig(cmtest_file)


#if __name__ == '__main__':
#    main()
