from bemio.io.aqwa import AqwaOutput
from bemio.io.output import write_hdf5

# Load AQWA output data file
aqwa_data = AqwaOutput(hydro_file='./data/aqwa_example_data.AH1',list_file='./data/aqwa_example_data.LIS')

# Calculate IRF and state space coefficients
for i in xrange(aqwa_data.data[0].num_bodies):
	aqwa_data.data[i].calc_irf_radiation(t_end=50,n_t=101,n_w=101)
	aqwa_data.data[i].calc_irf_excitation()
	# aqwa_data.data[i].calc_ss_radiation(max_order=3, r2_thresh=0.5 )
	

# Write hydrodynamic data to HDF5 file format
write_hdf5(aqwa_data)