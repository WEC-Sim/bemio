#! /usr/bin/python

# bemio WAMIT tutorial

# Load the needed modules from bemio
from bemio.io.wamit import WamitOutput
from bemio.io.output import write_hdf5


# Load the data using the wamit module.
# After running this line,  a object wamit_data of <type 'WamitOutput'> is created containing
# the data read from the wec3.out file. The data for each of the three bodies of the floating
# three-body oscillating flap is loaded into wamit_data.data[0], wamit_data.data[1],
# and wamit_data.data[2]
wamit_data = WamitOutput(out_file='./wec3.out')


# Calculate IRF and plot using the wamit module
# The loop below loops through each of the bodies and calculates the impulse response function
# and state space coefficients. In this case, the impulse response function is calculated from
# 0 to 50 seconds by setting t_end=50, using 101 timesteps n_t=101, and using 201 frequency steps
# n_w = 201
for i in xrange(wamit_data.data[0].num_bodies): #wamit_data.data[0].num_bodies
 	wamit_data.data[i].calc_irf_radiation()
	#wamit_data.data[i].calc_ss_radiation(max_order=5, r2_thresh=0.90)
	wamit_data.data[i].calc_irf_excitation()


# Save the data in the hdf5 format.
write_hdf5(wamit_data)
