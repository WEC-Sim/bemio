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
w = wio.WamitOutput(out_file='./data/wec3.out')

# Calculate IRF and plot
for i in xrange(w.data[0].nBodies):
	w.data[i].calcIRF(t_end=100., n_t = 1001, n_w=1001)

# Save the data in HDF5 and pickle format
hd.writeHdf5(w.data,w.files['hdf5'])
hd.writePickle(w.data,w.files['pickle'])