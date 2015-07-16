#! /usr/bin/python

from bemio.io.nemoh import NemohOutput


# ### Load the bem data object from bemio.data_structures #
# This module provides the functionality to save the data in the SHDF

from bemio.io.output import write_hdf5

# ### Read the Nemoh simulation data
nemoh_data = NemohOutput(sim_dir='./data/two_body', cal_file='Nemoh.cal', results_dir='Results', mesh_dir='Mesh')
nemoh_data.read_hydrostatics(body_num=1,file='./data/two_body/Mesh/Hydrostatics_1.dat')
nemoh_data.read_kh(body_num=1, file='./data/two_body/Mesh/KH_1.dat')

for i in xrange(nemoh_data.data[0].num_bodies):
	nemoh_data.data[i].calc_irf_radiation(t_end=50, n_t=101, n_w=201)
	nemoh_data.data[i].calc_irf_excitation(n_t=101, n_w=201)
	# nemoh_data.data[i].calc_ss_radiation()

write_hdf5(nemoh_data)
