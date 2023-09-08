#! /usr/bin/python

from bemio.io.nemoh import read

from bemio.io.output import write_hdf5

### Read the Nemoh simulation data
# One body
nemoh_data_one_body = read(sim_dir='./data/one_body', cal_file='Nemoh.cal', results_dir='Results', mesh_dir='Mesh')
nemoh_data_one_body.read_hydrostatics(body_num=0,file='./data/one_body/Mesh/Hydrostatics.dat')
nemoh_data_one_body.read_kh(body_num=0, file='./data/one_body/Mesh/KH.dat')
write_hdf5(nemoh_data_one_body)


# Two body
nemoh_data_two_body = read(sim_dir='./data/two_body', cal_file='Nemoh.cal', results_dir='Results', mesh_dir='Mesh')
nemoh_data_two_body.read_hydrostatics(body_num=1,file='./data/two_body/Mesh/Hydrostatics_1.dat')
nemoh_data_two_body.read_kh(body_num=1, file='./data/two_body/Mesh/KH_1.dat')

for i in range(nemoh_data_two_body.body[0].num_bodies):
	nemoh_data_two_body.body[i].calc_irf_radiation(t_end=20., n_t=1001, n_w=501)
	nemoh_data_two_body.body[i].calc_irf_excitation(t_end=20., n_t=1001, n_w=501)
	# nemoh_data_two_body.body[i].calc_ss_radiation()

write_hdf5(nemoh_data_two_body)
