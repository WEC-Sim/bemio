# Load the needed modules from bemio
from bemio.io.wamit import read
from bemio.io.output import write_hdf5


# Load the data using the wamit module.
wamit_data = read(out_file='wamit_data/oswec.out')


# Calculate IRF and plot using the wamit module
for i in range(wamit_data.body[0].num_bodies): #wamit_data.body[0].num_bodies
 	wamit_data.body[i].calc_irf_radiation(t_end=20, n_t=501, n_w=501, w_max=10.)
	#wamit_data.body[i].calc_ss_radiation(max_order=5, r2_thresh=0.90)
	wamit_data.body[i].calc_irf_excitation(t_end=30. , n_t=501, n_w=501)


# Save the data in the hdf5 format.
write_hdf5(wamit_data)
