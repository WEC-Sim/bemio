.. _getting_started:

Getting started
***************

.. _downloading:

Installing Python and the Required Dependencies
================================================
For Python Beginners
----------------------
It can be a pain to install Python, NumPy, SciPy, Matplotlib, h5py and other dependencies that are needed to run `bemio`. If you're new to Python, the easiest approach is to start by installing one of the science and engineering oriented Python distributions:

* `Anaconda <http://continuum.io/downloads>`_ (Linux, Mac, Windows)
* `PythonXY <https://code.google.com/p/pythonxy/>`_ (Windows)

.. Note::

	If you are using Anaconda, PythonXY or another non-Python.org distribution, you may need to install one of the packages identified in the `For Experienced Python Users`_ section. If you need to do so, you should instead install any needed modules via your distribution's own package manager.


For Experienced Python Users
-----------------------------------------
Install `Python 2.7.x <https://www.python.org/downloads/>`_ and the following Python packages using pip or your favorite package instillation method:
	
	* matplotlib
	* numpy
	* scipy
	* progressbar2
	* astropy
	* h5py

Installing HDFView
=======================
Although not required, `bemio` users will benefit from installing `HDFVIEW <http://www.hdfgroup.org/products/java/hdfview/>`_. _`HDFVIEW`  allows `bemio` users to view and interact with the standardized bemio hydrodynamic data files.

Downloading bemio
====================
`bemio` is distbuited through the `bemio GitHub web page <https://github.com/WEC-Sim/bemio/>`_. There are three ways to obtain the code:

Clone with Git (Recommended `bemio` for Users)
------------------------------------------------
`bemio` can obtained by cloning the repository with Git::

	git clone https://github.com/WEC-Sim/bemio.git

This method is recommended for most users because it makes it easy to update your local version of `bemio` to the latest version using Git's ``pull`` command::

	git pull

Fork with Git (Recommended for `bemio` Developers)
------------------------------------------------------------
If you are planning to contribute to the `bemio` code base, please follow the `forking instructions provided by GitHub <https://help.github.com/articles/fork-a-repo/>`_. If you implement an improvement you would like included in the `bemio` code base please make a `pull request <https://help.github.com/articles/using-pull-requests/>`_ so that your improvement can be merged into the code base.

Download code archive
------------------------
The easiest way to obtain a copy of `bemio` is to download the code archive:

* `bemio-1.0 <>`_. FIX THIS LINK.

If you chose this method you will have to re-download the code in order to receive code updates.

Installing `bemio`
==================
`bemio` is written in Python 2.7 and is compatible with OSX, Windows, and Linux operating systems. This section provides instructions on how to install `bemio` and test that you installation is working correctly. Simple instillation instructions are provided for users who simply want to run `bemio` and developer Instillation Instructions are provided for users who would like to develop and contribute to the `bemio` code base.

.. Note::
	
	In this section the location of the `bemio` source folder is refereed to as ``$BEMIO_SOURCE``

.. _simple_install:

Simple Installation Instructions 
--------------------------------

Step 1: Install Python 2.7.x
	Install Python 2.7.x as described in the `For Python Beginners`_ section

Step 2: Download `bemio`
	Download `bemio` as described in the `Downloading bemio`_ section

Step 3: Install `bemio`
	Open a command window (Windows) or shell (OSX and LInux) and navigate to ``$BEMIO_SOURCE``. Once inside the ``$BEMIO_SOURCE`` folder execute the following command::

		python setup.py install

	This installs the `bemio` files to your `$PYTHONPATH <https://docs.python.org/2/using/cmdline.html#environment-variables>`_. Note the above command must be run with administrative privileges on most systems.

Step 4: Test the instillation
	Test the `bemio` by running the following command from a command prompt inside the ``$BEMIO_SOURCE`` folder::

		python test.py

	If the  WAMIT, AQWA, and NEMOH tutorial cases (located in ``$BEMIO_SOURCE/tutorials``)  run successfully, you will receive get a success messages for each case.


.. _developer_install:

Developer Instillation Instructions
-----------------------------------

These instructions assume the user has a advanced level of Python knowledge. These instillation instructions should also be followed for users who do not have administrative rights on their system.

Step 1: Install Python 2.7.x
	Install Python 2.7.x and the Python modules identified in the `For Experienced Python Users`_ section.

Step 2: Download `bemio`
	Clone or fork `bemio` as described in the `Downloading bemio`_ section.

Step 3: Install `bemio`
	Add the ``$BEMIO_SOURCE`` folder to your `PYTHONPATH <https://docs.python.org/2/using/cmdline.html#environment-variables>`_ environment variable.

Step 4: Test the instillation
	Test the `bemio` by running the following command from a command prompt inside the ``$BEMIO_SOURCE`` folder::

		python test.py

	If the  WAMIT, AQWA, and NEMOH tutorial cases (located in ``$BEMIO_SOURCE/tutorials``)  run successfully, you will receive get a success messages for each case.
