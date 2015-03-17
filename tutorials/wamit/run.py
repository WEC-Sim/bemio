"""
Created on Mon Nov 17 23:51:21 2014

@author: mlawson

This is an example of how to use the wamitio module
"""

from bemio.io import wamit as wio
from bemio.data import bem as hd

# Load the data
w = wio.WamitOutput(out_file='./data/rm3_new.out')

# Calculate IRF and plot
for i in xrange(w.data[0].nBodies):
	w.data[i].calcIRF()


# Save the data in HDF5 and pickle format
hd.writeHdf5(w.data,w.files['hdf5'])