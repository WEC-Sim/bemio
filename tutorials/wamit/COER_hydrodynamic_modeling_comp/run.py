# Load the needed modules from bemio
from bemio.io.wamit import read
from bemio.io.output import write_hdf5

# Load the data using the wamit module.
wamit_data_f = read(out_file='wamit_data/coer_comp_f.out')

# Calculate IRF and SS coefficients
for i in range(wamit_data_f.body[0].num_bodies):
 	wamit_data_f.body[i].calc_irf_radiation(t_end=20.0, n_t=1001, n_w=501)
	wamit_data_f.body[i].calc_irf_excitation(t_end=20.0, n_t=1001, n_w=501)

# Save the data in the hdf5 format.
write_hdf5(wamit_data_f,out_file='coer_comp_f.h5')


# Load the data using the wamit module, scaling the data by rho=1000, g=9.81 and
# frequench (w)
wamit_data_f_scaled = read(out_file='wamit_data/coer_comp_f.out', scale=True)

# Calculate IRF and plot using the wamit module
for i in range(wamit_data_f_scaled.body[0].num_bodies):
 	wamit_data_f_scaled.body[i].calc_irf_radiation(t_end=20.0, n_t=1001, n_w=1001)
	wamit_data_f_scaled.body[i].calc_irf_excitation(t_end=20.0, n_t=1001, n_w=1001)

# Save the data in the hdf5 format.
write_hdf5(wamit_data_f_scaled,out_file='coer_comp_f_scaled.h5')
