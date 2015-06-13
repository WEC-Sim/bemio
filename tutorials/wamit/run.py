#! /usr/bin/python

# bemio WAMIT tutorial

# Import matplotlib - non required
import matplotlib.pyplot as plt
plt.interactive(True)

# Load the wamit module from bemio.io
# This module provides functionality to read, dimensionalize, visualizes the data.
from bemio.io import wamit


# Load the bem data object from bemio.data_structures
# This module provides the functionality to save the data in the SHDF
from bemio.io.hdf5 import write_hdf5


# Load the data using the wamit module.
# After running this line,  a object wamit_data of <type 'WamitOutput'> is created containing
# the data read from the wec3.out file. The data for each of the three bodies of the floating 
# three-body oscillating flap is loaded into wamit_data.data[0], wamit_data.data[1],
# and wamit_data.data[2]
wamit_data = wamit.WamitOutput(out_file='./data/rm3/rm3.out')


# Calculate IRF and plot using the wamit module
# The loop below loops through each of the bodies and calculates the impulse response function
# and state space coefficients. In this case, the impulse response function is calculated from
# 0 to 50 seconds by setting t_end=50, using 101 timesteps n_t=101, and using 201 frequency steps
# n_w = 201
for i in xrange(wamit_data.data[0].num_bodies):
	wamit_data.data[i].calc_irf_radiation(t_end=50,n_t=101,n_w=101)
	wamit_data.data[i].calc_ss_radiation()
	wamit_data.data[i].calc_irf_excitation(t_end=50,n_t=101,n_w=101)

# Save the data in SHDF format using the bem module
# The data is written to the file name specified by the out_file string. The wamit_data
# contains a dictionary named 'files' with an entry 'hdf5' which contains the default file name
# wec3.h5
write_hdf5(wamit_data)