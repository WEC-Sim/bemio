
# coding: utf-8

# ## NEMOH
# The following code describes the steps involved in using ``bemio`` to read, process, and save NEMOH coefficients. This code can be viewed in ``$BEMIO_SOURCE/tutorials/nemon/run.py``. An ipython notebook with the code is also provided in ``$BEMIO_SOURCE/tutorials/nemon/run.ipynb``.


# ### Load the nemoh module from bemio.io #
# This module provides functionality to read, dimensionalize, visualizes the data.

# In[ ]:

from bemio.io import nemoh


# ### Load the bem data object from bemio.data_structures #
# This module provides the functionality to save the data in the SHDF

# In[ ]:

from bemio.data_structures import bem


# ### Read the Nemoh simulation data from the 

# In[ ]:

nemoh_data = nemoh.NemohOutput(sim_dir='./data/two_body', cal_file='Nemoh.cal', results_dir='Results', mesh_dir='Mesh', out_name='./data/two_body.out')


# ### Read the hydrostatic and IH files
# These files must be read individually for each body as shown below

# In[ ]:

nemoh_data.read_hydrostatics(body_num=0,file='./data/two_body/Mesh/Hydrostatics_0.dat')
nemoh_data.read_kh(body_num=0, file='./data/two_body/Mesh/KH_0.dat')
nemoh_data.read_hydrostatics(body_num=1,file='./data/two_body/Mesh/Hydrostatics_1.dat')
nemoh_data.read_kh(body_num=1, file='./data/two_body/Mesh/KH_1.dat')


# ### Calculate the IRF and state space coefficients
# t_end, n_t and n_w are set to small numbers so the calculations run quickly. It is recommended that numbers are changed to 100, 1001, and 1001, respectively, for best accuracy of the IRF and state space coefficients.

# In[ ]:

for i in xrange(nemoh_data.data[0].num_bodies):
	nemoh_data.data[i].calc_irf(t_end=50, n_t=101, n_w=201)
	nemoh_data.data[i].calc_ss()


# ### Plot the surge and heave components of added mass and radiation damping for body 0

# In[ ]:

comps_to_plot = [[0,0],[2,2]]
nemoh_data.data[0].plot_am_rd(comps_to_plot)


# ### Write the data to the bemio data file format for use with WEC-Sim

# In[ ]:

bem.write_hdf5(nemoh_data)

