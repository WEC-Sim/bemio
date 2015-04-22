#! /usr/bin/python

# bemio AQWA tutorial

# Import matplotlib - non required
import matplotlib.pyplot as plt

# Load the aqwa module from bemio.io #
# This module provides functionality to read, dimensionalize, visualizes the data.
from bemio.io import aqwa

# Load the bem data object from bemio.data_structures #
# This module provides the functionality to save the data in the SHDF
from bemio.data_structures import bem

# Load AQWA output data file
aqwa_data = aqwa.AqwaOutput(out_file='./data/aqwa_example_data.lis')

# Plot diag components of added mass and damping
comps_to_plot = [[0,0],[1,1],[2,2]]
aqwa_data.data[0].plot_am_rd(comps_to_plot)
aqwa_data.data[1].plot_am_rd(comps_to_plot)
aqwa_data.data[2].plot_am_rd(comps_to_plot)

# Calculate IRF and plot using the aqwa module
# The loop below loops through each of the bodies and calculates the impulse response function
# and state space coefficients. In this case, the impulse response function is calculated from
# 0 to 50 seconds by setting t_end=50, using 101 timesteps n_t=101, and using 201 frequency steps
# n_w = 201
for i in xrange(aqwa_data.data[0].num_bodies):
	aqwa_data.data[i].calc_irf(t_end=50, n_t = 101, n_w=201)
	aqwa_data.data[i].calc_ss()

# Write hydrodynamic data to HDF5 file format
bem.write_hdf5(aqwa_data)
bem.write_pickle(aqwa_data)

# Close all plots
plt.close('all')