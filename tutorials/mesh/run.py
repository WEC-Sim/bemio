# Import the bemio.mesh_utilities module
from bemio.mesh_utilities import mesh
import numpy as np

# Read mesh
# non_symmetrical = mesh.read(file_name='./data/non_symmetrical.gdf')
# data = mesh.read(file_name='./data/buoy_a.vtp')
data = mesh.read(file_name='./data/ellipsoid.GDF')
data.center_of_gravity = np.array([0, 0, 0])

# Write meshes - uncomment these lines to wrtie the mesh to a different format
# data.write(mesh_format='VTK')
# data.write(mesh_format='NEMOH')
# data.write(mesh_format='WAMIT')

# View mesh - uncomment this line to view the mesh
# data.view()
