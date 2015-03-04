WAMIT
-----
WAMIT tutorial files can be found in the ``$BEMIO_TUTORIALS/wamit`` folder. As shown below, the folder contains a ``data`` folder and a ``run.py`` folder as shown below.

.. figure::  _static/wamit_tut_dir.png
   :align: center
   :width: 400pt

The data folder contains several example WAMIT output files. In this tutorial we will use the ``wec3.out`` example data. This file contains WAMIT simulation from a floating three-body oscillating flap device shown below. This device geometry was used in the Wave Energy Converter Code Comparison (WEC3) project [1] and described by Babarit et al. [2] (FIX CITATIONs).

.. figure::  _static/wec3_device.png
   :align: center
   :width: 400pt

The ``run.py`` file contains a custom written script that uses the `bemio` :ref:`API` to read the ``wec3.out`` file, plot hydrodynamic coefficients, and save the data in the SHDF.

.. literalinclude:: ../../tutorials/wamit/run.py

Running the script above produces the following output:

.. figure:: _static/wamit_tut_run.png
   :align: center
   :width: 500pt

.. figure:: _static/wamit_tut_fig.png
   :align: center
   :width: 400pt

Viewing the SHDF File
~~~~~~~~~~~~~~~~~~~~~

The SHDF file containing the WAMIT data was written to ``$BEMIO_TUTORIALS/wamit/data/wec3.h5`` and can be viewed using `HDFVIEW <http://www.hdfgroup.org/products/java/hdfview/>`_ as shown below. `HDFVIEW <http://www.hdfgroup.org/products/java/hdfview/>`_ allows the different data to be viewed and simple line plots to be made. For example, the image below shows the added mass matrix (``wamit_data_obj.data[0].am.all``) with a shape of 6 x 6*num_bodies x num_frequencies, the body name from the WAMTI file (``wamit_data_obj.data[0].name``), the body center of gravity from the WAMIT file (``wamit_data_obj.data[0].cg``). Also shown is the 1,1 (surge) component of the added mass matrix, which is stored in the SHDF file in the ``bodyN/added_mass/comps/comp_1_1`` file. `HDFVIEW <http://www.hdfgroup.org/products/java/hdfview/>`_ was used to plot this added mass component.

.. figure:: _static/wamit_tut_shdf.png
   :align: center
   :width: 550pt