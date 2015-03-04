#wamitio (WAMIT Input/Output)

##Metadata
* Author: Michael Lawson
* Date updated: 23 Feb 2015

##This module provides the functionality to
* Read hydrodynamic coefficients from WAMIT ".out" files
* Plot hydrodynamic coefficients
* Write hydrodynamic data to .mat files for use in WEC-SIm
* Write hydrodynamic data to HDF5 file format for use in other applications

##Description of files in this folder
* wamit.py: python module
* example/wamit-example.py: example of how to use the awqaio module
* example/oswec.out: example data of an WAMIT .out file from a simulation with three bodies
* example/skewed-sphere.out: example data of an WAMIT .out file from a simulation with one

##Notes
* This module is currently under active development, bug identification and fixed from the community are weclcome! Please post bugs and feature requests on GitHub.

##Requirements
* python 2.7 with the following packages:
  * numpy
  * scipy
  * matplotlib
  * pickle
  * h5py
