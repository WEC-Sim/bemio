from bemio.io import aqwa as aio
from bemio.io.output import write_hdf5

# Load AQWA output data file
aq = aio.AqwaOutput(hydro_file='./data/aqwa_example_data.AH1',list_file='./data/aqwa_example_data.LIS')

# Calculate IRF and state space coefficients
for i in xrange(aq.data[0].num_bodies):
	aq.data[i].calc_irf_radiation(t_end=50,n_t=101,n_w=101)
	aq.data[i].calc_ss_radiation(max_order=3, r2_thresh=0.5 )
	aq.data[i].calc_irf_excitation()

# Write hydrodynamic data to HDF5 file format
write_hdf5(aq)