.. _bemio_doc:

`bemio` Documentation and Users Guide
=====================================

.. _introduction:

Introduction
------------
The hydrodynamic coefficients (or data) that describe the radiation damping, added mass, wave diffraction force, and wave excitation force on floating bodies are widely used in multi body dynamics simulations of wave energy converters, ships, offshore platforms, and other floating structures to model hydrodynamic forces. These hydrodynamics coefficients are typically determined using frequency domain boundary element method (BEM) hydrodynamics codes such as `NEMOH <http://lheea.ec-nantes.fr/doku.php/emo/nemoh/start>`_, `WAMIT <http://www.wamit.com/>`_, and `AQWA <http://www.ansys.com/Products/Other+Products/ANSYS+AQWA>`_. `bemio` is a pre- and post-processing tool that works with these three codes. to reads the output data from these codes. `bemio` provides the functionality to:

* Read, view, and convert between different BEM code mesh formats. Specifically, `bemio` can read STL, VTK, WAMIT, and NEMOH mesh formats and convert between these formats.
* Read NEMOH, WAMIT, and AQWA simulation output files and save the data in a standardized human readable format that uses the `Hierarchical Data Format 5 (HDF5) <http://www.hdfgroup.org/HDF5/>`_.
* Calculate the impulse response functions (IRFs) and state space realization coefficients that represent the IRFs.

.. _motivation_for_devel:

Motivation for Development
--------------------------
There were two main motivations for developing `bemio`:

#. `bemio` was developed to support the WEC-Sim (Wave Energy Converter Simulator) project. During development and use of WEC-Sim, it became apparent WEC-Sim needed the ability to use hydrodynamic data from any BEM code and to allow users to input custom hydrodynamic data. Developing a standard hydrodynamic data input format for the code makes it easy for BEM data readers for different codes to be written and also provides a simple way for users to modify or specify custom hydrodynamic data.
#. Is is difficult to analyze, visualize, and compare hydrodynamic data, especially if the data is from different sources, Creating a standard significantly decreases the difficulty of these tasks.

.. _developers:
Developers
----------
* Michael Lawson, National Renewable Energy Laboratory
* Carlos Michelen, Sandia National Laboratories
* Yi-Hsiang Yu, National Renewable Energy Laboratory

Contents
--------
.. toctree::
   :maxdepth: 2
   
   license.rst
   installing.rst
   vtk.rst
   api_io.rst
   api_mesh.rst
   tutorials.rst
   module_doc.rst


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`