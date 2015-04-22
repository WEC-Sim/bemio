

from bemio.io import aqwa as aio
from bemio.data_structures import bem
import matplotlib.pyplot as plt
plt.close('all')
plt.interactive(True)

# Load AQWA output data file
aq = aio.AqwaOutput(out_file='./data/aqwa_example_data.AH1')

# Plot diag components of added mass and damping
comps_to_plot = [[0,0],[1,1],[2,2]]
aq.data[0].plot_am_rd(comps_to_plot)
aq.data[1].plot_am_rd(comps_to_plot)
aq.data[2].plot_am_rd(comps_to_plot)

# Calculate IRF and plot
for i in xrange(aq.data[0].num_bodies):
	aq.data[i].calc_irf()
	#aq.data[i].calc_ss()

# Write hydrodynamic data to HDF5 file format
bem.write_hdf5(aq.data,aq.files['hdf5'])
bem.write_pickle(aq.data,aq.files['pickle'])

# Keep Python running so the user can view the plots #
raw_input('\nPress enter to exit and close plots')
plt.close('all')

