Using the bemio.io Module
=========================
This section describes how to build custom data readers for NEMOH, WAMIT, and AQWA output data using the bemio.io module.

Loading Data
~~~~~~~~~~~~
Users can load data by interacting with the following with the following bemio modules

	* :mod:`bemio.io.nemoh`
	* :mod:`bemio.io.wamit`
	* :mod:`bemio.io.aqwa`

The Python code below shows how a user would use bemio to read NEMOH, WAMIT, and AQWA data. Notice that the :mod:`bemio.io` objects are used to create data objects that contain the hydrodynamic data from the data files produced by the hydrodynamic simulation codes.::

	from bemio.io.nemoh import read as read_nemoh  # Import NEMOH module
	from bemio.io.wamit import read as read_wamit  # Import WAMIT module
	from bemio.io.aqwa import read as read_aqwa   # Import AQWA module

	# Read  NEMOH data into the nemoh_data object
	nemoh_data_obj = read_nemoh(sim_dir='$BEMIO_SOURCE/tutorials/nemoh/data/')

	# Read data from a WAMIT .out file into the wamit_data object
	wamit_data_obj = read_wamit(out_file='$BEMIO_SOURCE/tutorials/wamit/rm3/wamit_data/rm3.out')

	# Read data from a AQWA .lis file into the aqwa_data object
	aqwa_data_obj = read_aqwa(hydro_file='$BEMIO_SOURCE/tutorials/aqwa/aqwa_example_data.AH1', list_file=''$BEMIO_SOURCE/tutorials/aqwa/aqwa_example_data.LIS')

.. _irf_and_ss:

Data Objects
~~~~~~~~~~~~
The bemio data objects described in the previous contain BEM simulation data stored in a standard format. Details on the bem data forma can be found here :mod:`bemio.data_structures.bem`

Calculating Impulse Response Functions and Sate Space Coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

	# Calculate the radiation damping impulse response functions
	nemoh_data_obj.data[i].calc_irf_radiation()

	# Calculate the state space coefficients corresponding to the radiation damping IRFs
	nemoh_data_obj.data[i].calc_ss_radiation()

	# Calculate the excitation force impulse response function
	nemoh_data_obj.data[i].calc_ss_excitation()


Writing the Data to the Standard bemio Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bemio provides a consistent format to store BEM simulation data from NEMOH, WAMT, AQWA. The bemio format uses the `Hierarchical Data Format 5 (HDF5) <http://www.hdfgroup.org/HDF5/>`_, a data structure that is similar to a computers file system, to store complex data in a simple format. HDF5 files can be easily viewed and modified with `HDF5VIEW <http://www.hdfgroup.org/products/java/hdfview/>`_.

`WEC-Sim <https://github.com/WEC-Sim/WEC-Sim>`_ uses the bemio format as an input file for hydrodynamic data used in wave energy converter simulations.

The bemio :meth:`bemio.io.output.write_hdf5` function can be used to write the data to a SHDF file::

	# Import the bem module
	from bemio.io.output import write_hdf5

	# Write the data to a SHDF file named shdf_example.h5
	write_hdf5(data=nemoh_data_obj.data)

Once the data is written to the HDF5 file, it can be viewed using `HDFVIEW`_. The figure below shows the structure of the SHDF produced by bemio from the NEMOH data. Note that the NEMOH simulation was of a two-body point absorber.

.. figure::  _static/rm3_hdf5.png
   :align:   center
   :width: 600pt

As shown in the figure above above, the SHDF contains three top levels:

	``bem_data``:
		This folder contains the raw and unprocessed output data files from the BEM code
	``bodyN``:
		The folders named ``body1``, ``body2``, ... ``bodyN`` contain the processed data for that body. In the case above there are two bodies, corresponding to the two bodies of the two body point absorber. Notice that the folders beneath each ``bodyN`` folder contains the hydrodynamic data for that body.
	``simulatino_parameters``:
		This data structure contains simulation parameters from the BEM simulation that are independent of the body number. For example, the ``water_depth`` variable contains the water depth used during the BEM simulation
