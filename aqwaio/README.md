#Metadata
* author: michael lawson
* date updated: 6 Feb 2015

#This module provides the functionality to
* Read hydrodynamic coefficients from AQWA ".lis" files
* Plot hydrodynamic coefficients
* Write hydrodynamic data to .mat files for use in WEC-SIm
* Write hydrodynamic data to HDF5 file format for use in other applications

#Description of files in this folder
* aqwaio.py: python module
* example/aqwaio-example.py: example of how to use the awqaio module
* example/aqwa-example-data.lis: example data of an aqwa run. no hydrodynamic interactions for a three body flapping device that consists of a frame and two flaps
* example/compareBEM.m: a matlab file used during development - normal users can ignore this file
aqwa-manually-generated.mad: a matlab file used during development - normal users can ignore this file

#Notes
* Only works for AQWA simulations without body to body interactions turned on
* This module is currently under active development, bug identification and fixed from the community are weclcome! Please post bugs and feature requests on GitHub.

#Requirements
* python 2.7 with the following packages:
  * numpy
  * scipy
  * matplotlib
  * h5py
