# Import modules
import numpy as np
import h5py
from bemio.data_structures.wave_excitation import convolution

# Load the fast data
fast = np.loadtxt('OMAERB.out',skiprows=8)

# Load wave time series
wave_data = np.loadtxt('wave_elevation.dat')
eta = wave_data[:,1]/1000. # Convert from mm to m
eta_t = wave_data[:,0]

# Load excitation force IRF
hydro_data_f = h5py.File('../../wamit/COER_hydrodynamic_modeling_comp/coer_comp_f_scaled.h5','r+')
irf_f = np.array(hydro_data_f['/body1/hydro_coeffs/excitation/impulse_response_fun/f'])[2,0,:]
irf_t_f = np.array(hydro_data_f['/body1/hydro_coeffs/excitation/impulse_response_fun/t'])

ex_f = convolution(irf=irf_f, irf_t=irf_t_f, eta=eta, eta_t=eta_t)

# Plot the data against the results from the NREL FAST code for verification.
# Set interactive and close all active plots
try:
    import matplotlib.pyplot as plt
    plt.interactive(True)
    plt.close('all')
    plt.figure()
    plt.plot(ex_f.excitation_force.t, ex_f.excitation_force.f, 'r--',label='bemio')
    plt.plot(fast[:,0],fast[:,56],label='FAST')
    plt.legend()
except:
    pass
