
from bemio.io import aqwa as aio
from bemio.data_structures import bem
import matplotlib.pyplot as plt

# Load AQWA output data file
aq = aio.AqwaOutput(out_file='./data/aqwa_example_data.AH1')

# Plot diag components of added mass and damping
comps_to_plot = [[0,0],[1,1],[2,2]]
aq.data[0].plot_am_rd(comps_to_plot)
aq.data[1].plot_am_rd(comps_to_plot)
aq.data[2].plot_am_rd(comps_to_plot)

# Calculate IRF and plot using the aqwa module
#for i in xrange(aq.data[0].num_bodies):
	#aq.data[i].calc_irf(t_end=50, n_t = 101, n_w=201)
	#aq.data[i].calc_ss()

# Write hydrodynamic data to HDF5 file format
bem.write_hdf5(aq)
bem.write_pickle(aq)

# Keep Python running so the user can view the plots #
plt.show()
raw_input('\nPress enter to exit and close plots')
plt.close('all')

