# Import modules
import numpy as np
import h5py
#import matplotlib.pyplot as plt
from bemio.data_structures.wave_excitation import convolution

# Set interactive and close all active plots
#plt.interactive(True)
#plt.close('all')

# Load the fast data
fast = np.loadtxt('OMAERB.out',skiprows=8)

# Load wave time series
wave_data = np.loadtxt('wave_elevation.dat')
eta = wave_data[:,1]/1000. # Convert from mm to m
eta_t = wave_data[:,0]

# Load excitation force IRF
hydro_data_f = h5py.File('coer_comp_scaled_f.h5','r+')
irf_f = np.array(hydro_data_f['/body1/hydro_coeffs/excitation/impulse_response_fun/f'])[2,0,:]
irf_t_f = np.array(hydro_data_f['/body1/hydro_coeffs/excitation/impulse_response_fun/t'])

# hydro_data_f2 = h5py.File('coer_comp_scaled_f2.h5','r+')
# irf_f2 = np.array(hydro_data_f2['/body1/hydro_coeffs/excitation/impulse_response_fun/f'])[2,0,:]
# irf_t_f2 = np.array(hydro_data_f2['/body1/hydro_coeffs/excitation/impulse_response_fun/t'])

# Calculate excitation force coefficients and the wave elevation using the
# convolution method to deterine the excitation force time series
ex_f = convolution(irf=irf_f, irf_t=irf_t_f, eta=eta, eta_t=eta_t)
# ex_f2 = convolution(irf=irf_f2, irf_t=irf_t_f2, eta=eta, eta_t=eta_t)

# Plot the data against the results from the NREL FAST code for verification.
#plt.figure()
#plt.plot(ex_f.excitation_force.t,ex_f.excitation_force.f,'r--',label='bemio_conv_f')
# plt.plot(ex_f2.excitation_force.t,ex_f2.excitation_force.f,'g--',label='bemio_conv_f2')
#plt.plot(fast[:,0],fast[:,56],label='FAST')
#plt.legend()
