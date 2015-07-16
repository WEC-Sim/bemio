# Import the bemio.mesh_utilities module
from bemio.mesh_utilities import mesh
import numpy as np
from copy import copy

# Read mesh
ellipsoid = mesh.read(file_name='ellipsoid.gdf')

# Creat a copy of the mesh
ellipsoid_scale_translate = copy(ellipsoid)

# Scale and translate the copy of the mesh
ellipsoid_scale_translate.scale(scale=[0,10.,.5])
ellipsoid_scale_translate.translate(translate=[1.,2.,3.])

# Write the scaled and tranalated mesh to NEMOH format
ellipsoid_scale_translate.write(mesh_format='NEMOH')
