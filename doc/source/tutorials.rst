.. toctree::

Tutorials
=========
This section provides tutorials that show how to use :ref:`API` to read, process, and save NEMOH, WAMIT, and AQWA simulation data. Specifically, the tutorials show how to:

	#. Create custom programs that read BEM simulation data
	#. Calculate impulse response functions
	#. Calculate state space coefficients from the impulse response functions
	#. Plot the data for visualization
	#. Save the data to the SHDF

`bemio` was developed so that users interact with the API in the same way regardless of the type of BEM data that is being read and processed. As such, we will present a thorough description of how to process WAMIT data in the 'WAMIT'_ tutorial, and abbreviated versions of the NEMOH and AQWA tutorials.


.. _data_files:

Tutorial Data Files
-------------------
Tutorial files are distributed with the `bemio` code and can be found in the ``$BEMIO_SOURCE/tutorials`` folder. This location will be refereed to as ``$BEMIO_TUTORIALS``.

.. _wamit_tut:


Available Tutorials
------------------------
.. toctree::
   :maxdepth: 1
   
   wamit.rst
   nemoh.rst
   aqwa.rst
