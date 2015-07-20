#! /usr/bin/python

# Import the bemio.mesh_utilities module
from bemio.mesh_utilities import mesh
import numpy as np
from copy import copy

# Read mesh
ellipsoid = mesh.read(file_name='ellipsoid.gdf')

# Scale and translate the copy of the mesh
ellipsoid.scale(scale_vect=[5.,10.,10.])
ellipsoid.translate(translation_vect=[1.,2.,3.])

# Write the scaled and tranalated mesh to WAMIT format
ellipsoid.write(mesh_format='WAMIT')
