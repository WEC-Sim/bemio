# Import the bemio.mesh_utilities module
from bemio.mesh_utilities import mesh
import numpy as np

# Read mesh
# non_symmetrical = mesh.read(file_name='./data/non_symmetrical.gdf')
data = mesh.read(file_name='./data/ellipsoid.GDF')
data.center_of_gravity = np.array([0, 0, 0])
# Write meshes
# data.write(mesh_format='VTK')
# data.write(mesh_format='NEMOH')
# data.write(mesh_format='WAMIT')

# View mesh
# data.view()
