"""
Created on Thu Jan 15 10:42:54 2015

@author: mlawson
"""
from bemio.io import aqwa as aio
from bemio.data import bem as hd
import matplotlib.pyplot as plt
plt.close('all')
plt.interactive(True)

# Load AQWA output data file
aq = aio.AqwaOutput(outFile='./data/aqwa_example_data.lis')

# Plot diag components of added mass and damping
comps_to_plot = [[0,0],[1,1],[2,2]]
aq.data[0].plot_am_rd(comps_to_plot)
aq.data[1].plot_am_rd(comps_to_plot)
aq.data[2].plot_am_rd(comps_to_plot)

# Write hydrodynamic data to HDF5 file format
hd.write_hdf5(aq.data,aq.files['hdf5'])
hd.write_pickle(aq.data,aq.files['pickle'])