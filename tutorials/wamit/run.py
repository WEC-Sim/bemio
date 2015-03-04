"""
Created on Mon Nov 17 23:51:21 2014

@author: mlawson

This is an example of how to use the wamitio module
"""

from bemio.io import wamit as wio
from bemio.data import bem as hd
import matplotlib.pyplot as plt
plt.close('all')
plt.interactive(True)

# Load the data
w = wio.WamitOutput(outFile='./data/oswec.out')

# Plot selected components of hydrodynamic coefficinets and excitation force
# Note that python uses zero indexing
comps = [[0,0],[1,1],[2,2]]
w.data[0].plotAddedMassAndDamping(comps)
w.data[1].plotAddedMassAndDamping(comps)
w.data[2].plotAddedMassAndDamping(comps)
w.data[0].plotExcitation([0])

# Save the data in HDF5 format
hd.writeHdf5(w.data,w.files['hdf5'])
hd.writePickle(w.data,w.files['pickle'])