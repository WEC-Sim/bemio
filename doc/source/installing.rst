Downloading and Installing
==========================

Downloading bemio
-----------------
bemio is distributed through the `bemio GitHub repositiory <https://github.com/WEC-Sim/bemio/>`_.

.. note::

	The folder that contains the bemio source code will be referred to as the ``$BEMIO_SOURE`` folder.

Dependencies
-------------
**Python:** `Python 2.7.x <https://www.python.org/downloads/>`_ and the following Python packages are required to run bemio. These packages can easily be installed using using `pip <https://pypi.python.org/pypi/pip>`_  or your preferred package instillation method:

	* matplotlib
	* numpy
	* scipy
	* progressbar2
	* astropy
	* h5py

**Visualization Toolkit (VTK):** If you would like to use the ``bemio.mesh_utilities`` you will also have to install VTK and the Python 2.7.x bindings. Instructions for how to install VTK are provided in the :ref:`Installing VTK with Python Bindings` section.

**HDFVIEW:** Although not required, bemio users will benefit from installing `HDFVIEW <http://www.hdfgroup.org/products/java/hdfview/>`_. HDFVIEW allows bemio users to view and interact with the bemio output files.

Installing: Quick Start for Python Beginners
--------------------------------------------
**Step 1 - Install Python:** It can be a pain to install Python, NumPy, SciPy, Matplotlib, h5py and other dependencies that are needed to run bemio. If you're new to Python, the easiest approach is to start by installing one of the science and engineering oriented Python distributions:

	* `Anaconda <http://continuum.io/downloads>`_ (Linux, Mac, Windows)
	* `PythonXY <https://code.google.com/p/pythonxy/>`_ (Windows)

.. Note::

	If you are using Anaconda, PythonXY or another non-Python.org distribution, you may need to install one of the packages identified in the `Dependencies`_ section. If you need to do so, you should install any needed modules via your distribution's package manager.

**Step 2 - Download bemio:** Download bemio from the `GitHub repositiory <https://github.com/WEC-Sim/bemio/>`_

**Step 3 - Install bemio:** Open a command window (Windows) or terminal window (OSX and Linux) and navigate to ``$BEMIO_SOURCE``. Once inside the ``$BEMIO_SOURCE`` folder execute the following command to install bemio::

	python setup.py install --user

**Step 4 - Test the instillation:** Test the bemio by running the following command from a command prompt inside the ``$BEMIO_SOURCE`` folder::

	python install_test.py

If the tutorial cases (located in ``$BEMIO_SOURCE/tutorials``) run successfully you will receive a success messages.


Installing: For Experienced Python Users
-----------------------------------------
These instructions assume the user has an advanced level of Python knowledge.

**Step 1 - Install Python 2.7.x:** Install Python 2.7.x and the Python modules identified in the `Dependencies`_ section.

**Step 2 - Download bemio:** Download bemio from the `GitHub repositiory <https://github.com/WEC-Sim/bemio/>`_

**Step 3 - Install bemio:** Add the ``$BEMIO_SOURCE`` folder to your `PYTHONPATH <https://docs.python.org/2/using/cmdline.html#environment-variables>`_.

**Step 4 - Test the instillation:** Test bemio by running the following command from a command prompt inside the ``$BEMIO_SOURCE`` folder::

	python install_test.py

If the tutorial cases (located in ``$BEMIO_SOURCE/tutorials``) run successfully you will receive a success messages.
