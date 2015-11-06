.. _bemio_doc:

`bemio` Documentation and Users Guide
=====================================

.. _introduction:

Introduction
------------
The hydrodynamic coefficients (or data) that describe the radiation damping, added mass, wave diffraction force, and wave excitation force on floating bodies are widely used in multi body dynamics simulations of wave energy converters, ships, offshore platforms, and other floating structures to model hydrodynamic forces. These hydrodynamics coefficients are typically determined using frequency domain boundary element method (BEM) hydrodynamics codes such as `NEMOH <http://lheea.ec-nantes.fr/doku.php/emo/nemoh/start>`_, `WAMIT <http://www.wamit.com/>`_, and `AQWA <http://www.ansys.com/Products/Other+Products/ANSYS+AQWA>`_. bemio is a pre- and post-processing tool that works with these three codes.

What can bemio do?
------------------

bemio provides the following functionality:
	* **Mesh I/O:** Read, view, and convert between mesh formats used by BEM codes. Specifically, bemio can read STL, VTK, WAMIT, and NEMOH mesh formats and convert between them.
	* **BEM code results I/O:** Read NEMOH, WAMIT, and AQWA simulation output files and save the data in a standardized human readable bemio format that uses the `Hierarchical Data Format 5 (HDF5) <http://www.hdfgroup.org/HDF5/>`_.
	* **Calculation of impulse response function** Calculate the wave excitation and radiation dampinf impulse response functions (IRFs).
	* **Calculation of state space realization coefficients:** Calculation of state space realization coefficients that represent the IRFs.

Developers
----------
* Michael Lawson, National Renewable Energy Laboratory
* Carlos Michelen, Sandia National Laboratories
* Yi-Hsiang Yu, National Renewable Energy Laboratory

Contents
--------
.. toctree::
   :maxdepth: 2

   installing.rst
   vtk.rst
   api_io.rst
   tutorials.rst
   module_doc.rst
   license.rst


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
