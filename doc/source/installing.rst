Downloading and Installing
==========================

Downloading bemio
-----------------
<<<<<<< HEAD
bemio is distbuited through the `bemio GitHub web page <https://github.com/WEC-Sim/bemio/>`_ and can be downloaded using one of the following methods.

**Clone with Git (recommended):** bemio can obtained by cloning the repository using the GitHub GUI (`Mac <https://mac.github.com/>`_ and `Windows <https://windows.github.com/>`_ ) or from the command line using the command::

 	git clone https://github.com/WEC-Sim/bemio.git

This method is recommended because it makes it easy to update your local version of bemio to the latest version using the ``git pull`` command.

**Download code archive:** The code archive can be downloaded from `bemio-1.0 <https://github.com/WEC-Sim/bemio/archive/master.zip>`_.
=======
bemio is distributed through the `bemio GitHub web page <https://github.com/WEC-Sim/bemio/>`_. The code can be obtained using one of the following methods:

	* **Clone with GitHub GUI:** bemio can obtained by cloning the repository using the GitHub GUI (available for `Mac <https://mac.github.com/>`_ and `Windows <https://windows.github.com/>`_ ).
	* **Clone from command line:** ``git clone https://github.com/WEC-Sim/bemio.git``
	* **Download code archive:** `bemio-1.x.x <https://github.com/WEC-Sim/bemio/archive/master.zip>`_.
>>>>>>> 66c9ac2ebca143fd02284d952959fe4dd806a834

.. note::
	
	The folder that contains the bemio source code will be referred to as the ``$BEMIO_SOURE`` folder.

Dependencies
-------------
<<<<<<< HEAD
Bemio requires the following software and packages.

**Python:** `Python 2.7.x <https://www.python.org/downloads/>`_ and the following Python packages are required to run WEC-Sim. These packages can easily be installed using using `pip <https://pypi.python.org/pypi/pip>`_  or your favorite package instillation method:
=======
**Python:** `Python 2.7.x <https://www.python.org/downloads/>`_ and the following Python packages are required to run bemio. These packages can easily be installed using using `pip <https://pypi.python.org/pypi/pip>`_  or your preferred package instillation method:
>>>>>>> 66c9ac2ebca143fd02284d952959fe4dd806a834

	* matplotlib
	* numpy
	* scipy
	* progressbar2
	* astropy
	* h5py

**Visualization Toolkit (VTK):** If you would like to use the ``bemio.mesh_utilities`` you will also have to install VTK and the Python 2.7.x bindings. Instructions for how to install VTK are provided in the :ref:`Installing VTK with Python Bindings` section.

<<<<<<< HEAD
**HDFVIEW:** Although not required, bemio users will benefit from installing `HDFVIEW <http://www.hdfgroup.org/products/java/hdfview/>`_. `HDFVIEW`  allows bemio users to view and interact with the bemio output files.
=======
**HDFVIEW:** Although not required, bemio users will benefit from installing `HDFVIEW <http://www.hdfgroup.org/products/java/hdfview/>`_. HDFVIEW allows bemio users to view and interact with the bemio output files.
>>>>>>> 66c9ac2ebca143fd02284d952959fe4dd806a834

Installing: Quick Start for Python Beginners 
--------------------------------------------
**Step 1 - Install Python:** It can be a pain to install Python, NumPy, SciPy, Matplotlib, h5py and other dependencies that are needed to run bemio. If you're new to Python, the easiest approach is to start by installing one of the science and engineering oriented Python distributions:
	
	* `Anaconda <http://continuum.io/downloads>`_ (Linux, Mac, Windows)
	* `PythonXY <https://code.google.com/p/pythonxy/>`_ (Windows)
	
.. Note::

	If you are using Anaconda, PythonXY or another non-Python.org distribution, you may need to install one of the packages identified in the `Dependencies`_ section. If you need to do so, you should install any needed modules via your distribution's package manager.

**Step 2 - Download bemio:** Download bemio as described in the `Downloading bemio`_ section

<<<<<<< HEAD
**Step 3 - Install bemio:** Open a command window (Windows) or shell (OSX and LInux) and navigate to ``$BEMIO_SOURCE``. Once inside the ``$BEMIO_SOURCE`` folder execute the following command::
=======
**Step 3 - Install bemio:** Open a command window (Windows) or terminal window (OSX and Linux) and navigate to ``$BEMIO_SOURCE``. Once inside the ``$BEMIO_SOURCE`` folder execute the following command to install bemio::
>>>>>>> 66c9ac2ebca143fd02284d952959fe4dd806a834

	python setup.py install --user

**Step 4 - Test the instillation:** Test the bemio by running the following command from a command prompt inside the ``$BEMIO_SOURCE`` folder::

	python test.py

If the  WAMIT, AQWA, and NEMOH tutorial cases (located in ``$BEMIO_SOURCE/tutorials``)  run successfully, you will receive get a success messages for each case.


Installing: For Experienced Python Users
-----------------------------------------
These instructions assume the user has an advanced level of Python knowledge.

**Step 1 - Install Python 2.7.x:** Install Python 2.7.x and the Python modules identified in the `Dependencies`_ section.

<<<<<<< HEAD
**Step 2 - Download bemio:** Clone or fork bemio from the `bemio GitHub page <https://github.com/WEC-Sim/bemio/>`_ using git::

	git clone https://github.com/WEC-Sim/bemio.git

**Step 3 - Install bemio:** Add the ``$BEMIO_SOURCE`` folder to your `PYTHONPATH <https://docs.python.org/2/using/cmdline.html#environment-variables>`_ environment variable.
=======
**Step 2 - Download bemio:** See the `Downloading bemio`_ section

**Step 3 - Install bemio:** Add the ``$BEMIO_SOURCE`` folder to your `PYTHONPATH <https://docs.python.org/2/using/cmdline.html#environment-variables>`_.
>>>>>>> 66c9ac2ebca143fd02284d952959fe4dd806a834

**Step 4 - Test the instillation:** Test the bemio by running the following command from a command prompt inside the ``$BEMIO_SOURCE`` folder::

	python test.py

If the  WAMIT, AQWA, and NEMOH tutorial cases (located in ``$BEMIO_SOURCE/tutorials``)  run successfully, you will receive get a success messages for each case.
