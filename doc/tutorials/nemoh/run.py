#! /usr/bin/python

# bemio NEMOH tutorial

# Import matplotlib and make the session interactive for easy plot visualization #
import matplotlib.pyplot as plt


# Load the nemoh module from bemio.io #
# This module provides functionality to read, dimensionalize, visualizes the data.
from bemio.io import nemoh


# Load the bem data object from bemio.data_structures #
# This module provides the functionality to save the data in the SHDF
from bemio.data_structures import bem

# Read the Nemoh simulation data from the 
nemoh_data = nemoh.NemohOutput(sim_dir='./data/two_body', cal_file='Nemoh.cal', results_dir='Results', mesh_dir='Mesh', out_name='./data/two_body.out')

# Read the hydrostatic and kH files
# These files must be read individually for each body as shown below
nemoh_data.read_hydrostatics(body_num=0,file='./data/two_body/Mesh/Hydrostatics_0.dat')
nemoh_data.read_kh(body_num=0, file='./data/two_body/Mesh/KH_0.dat')
nemoh_data.read_hydrostatics(body_num=1,file='./data/two_body/Mesh/Hydrostatics_1.dat')
nemoh_data.read_kh(body_num=1, file='./data/two_body/Mesh/KH_1.dat')

# Calculate the IRF and state space coefficients
# t_end, n_t and n_w are set to small numbers so the calculations run quickly. It is recommended that
# numbers are changed to 100, 1001, and 1001, respectively, for best accuracy of the IRF and state
# space coefficients.
for i in xrange(nemoh_data.data[0].num_bodies):
	nemoh_data.data[i].calc_irf(t_end=50, n_t=101, n_w=201)
	nemoh_data.data[i].calc_ss()

# Write the data to the bemio data file format for use with WEC-Sim
bem.write_hdf5(data=nemoh_data.data, out_file=nemoh_data.files['hdf5'])

show()

