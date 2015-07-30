import numpy as np
import h5py
import matplotlib.pyplot as plt
from bemio.data_structures.wave_excitation import convolution

# Set interactive and close all active plots
plt.interactive(True)
plt.close('all')

A = 0.025
w = (2.*np.pi)/3.
eta_t = np.linspace(0,200,1000)
eta = A*np.cos(w*eta_t)




# Load excitation force IRF
hydro_data_f = h5py.File('coer_comp_scaled_f.h5','r+')
irf_f = np.array(hydro_data_f['/body1/hydro_coeffs/excitation/impulse_response_fun/f'])[2,0,:]
irf_t_f = np.array(hydro_data_f['/body1/hydro_coeffs/excitation/impulse_response_fun/t'])

# Calculate excitation force using the convolution method
ex = convolution(irf=irf_f, irf_t=irf_t_f, eta=eta, eta_t=eta_t)

# Plot the data
plt.figure()
plt.plot(ex.excitation_force.t,ex.excitation_force.f,'r--',label='bemio_conv_f')
