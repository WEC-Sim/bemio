#! /usr/bin/python

# bemio WAMIT tutorial

# Import matplotlib - non required
import matplotlib.pyplot as plt

# Load the wamit module from bemio.io
# This module provides functionality to read, dimensionalize, visualizes the data.
from bemio.io import wamit


# Load the bem data object from bemio.data_structures
# This module provides the functionality to save the data in the SHDF
from bemio.data_structures import bem


# Load the data using the wamit module.
# After running this line,  a object wamit_data of <type 'WamitOutput'> is created containing
# the data read from the wec3.out file. The data for each of the three bodies of the floating 
# three-body oscillating flap is loaded into wamit_data.data[0], wamit_data.data[1],
# and wamit_data.data[2]
wamit_data = wamit.WamitOutput(out_file='./data/wec3.out')


# Plot the added hydrodynamic coefficients for visualization using the wamit module #
# Specifically, the surge (0,0 component) and heave (2,2 component) for body0 in the WAMTI
# simulation are plotted using matplotlib. Note that Python uses zero indexing.
components_to_plot = [[0,0],[2,2]]
wamit_data.data[0].plot_am_rd(components_to_plot)


# Calculate IRF and plot using the wamit module
# The loop below loops through each of the bodies and calculates the impulse response function
# and state space coefficients. In this case, the impulse response function is calculated from
# 0 to 50 seconds by setting t_end=50, using 101 timesteps n_t=101, and using 201 frequency steps
# n_w = 201
for i in xrange(wamit_data.data[0].num_bodies):
	wamit_data.data[i].calc_irf(t_end=50, n_t = 101, n_w=201)
	wamit_data.data[i].calc_ss()

# Save the data in SHDF format using the bem module
# The data is written to the file name specified by the out_file string. The wamit_data
# contains a dictionary named 'files' with an entry 'hdf5' which contains the default file name
# wec3.h5
bem.write_hdf5(wamit_data)


# Close all plots
plt.close('all')