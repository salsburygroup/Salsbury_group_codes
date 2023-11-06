#!/usr/env python
"""
An example script that uses distances between specified atoms as the clustering input. This usage also shows the need
to eventually add a featurizer module.
"""

from ClusterPrediction.DCDReader import DCDReader
from ClusterPrediction.Clusterer import HDBSCAN
from ClusterPrediction.Clusterer import IMWKRescaled
import ClusterPrediction.Plotter as Plotter
import ClusterPrediction.Saver as Saver
import mdtraj as md
import numpy as np

t = DCDReader(
    '/Volumes/RyanMdata/F10/ZnCl2_150mM/cat16/F10andZincOnly.dcd',
    '/Volumes/RyanMdata/F10/ZnCl2_150mM/cat16/F10andZincOnly.pdb').load()

F5_Zn = t.topology.select_pairs('element F', 'element Zn')
FZd = md.compute_distances(t, F5_Zn).astype('float64')

O4_Zn = t.topology.select_pairs('name O4', 'element Zn')
O4Zd = md.compute_distances(t, O4_Zn).astype('float64')

O2_Zn = t.topology.select_pairs('name O2', 'element Zn')
O2Zd = md.compute_distances(t, O2_Zn).astype('float64')

P_Zn = t.topology.select_pairs('element P', 'element Zn')
PZd = md.compute_distances(t, P_Zn).astype('float64')

features = np.concatenate((FZd, O4Zd, O2Zd, PZd), axis=1).astype('float64')

HDBSCAN_labels = HDBSCAN(features).fit()
Amorim_labels = IMWKRescaled(features).fit()

Saver.TimeSeries("/Volumes/RyanMdata/F10/ZnCl2_150mM/cat16/Analysis/Clusters/withIons/F5O4O2P/HDBSCAN/labels.txt",
                 HDBSCAN_labels).save()
# Plot time series
Plotter.TimeSeries("/Volumes/RyanMdata/F10/ZnCl2_150mM/cat16/Analysis/Clusters/withIons/F5O4O2P/HDBSCAN/timeseris.png",
                   HDBSCAN_labels).plot()

Saver.PDB("/Volumes/RyanMdata/F10/ZnCl2_150mM/cat16/Analysis/Clusters/withIons/F5O4O2P/HDBSCAN/",
          HDBSCAN_labels, t, atom_selection='not element H')

Saver.TimeSeries("/Volumes/RyanMdata/F10/ZnCl2_150mM/cat16/Analysis/Clusters/withIons/F5O4O2P/Amorim/labels.txt",
                 Amorim_labels).save()
# Plot time series
Plotter.TimeSeries("/Volumes/RyanMdata/F10/ZnCl2_150mM/cat16/Analysis/Clusters/withIons/F5O4O2P/Amorim/timeseris.png",
                   Amorim_labels).plot()

Saver.PDB("/Volumes/RyanMdata/F10/ZnCl2_150mM/cat16/Analysis/Clusters/withIons/F5O4O2P/Amorim/",
          Amorim_labels, t, atom_selection='not element H')
