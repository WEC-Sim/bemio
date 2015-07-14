# Import the bemio.mesh_utilities module
from bemio.mesh_utilities import mesh
import numpy as np

# Read mesh
# non_symmetrical = mesh.read(file_name='./data/non_symmetrical.gdf')
# data = mesh.read(file_name='./data/buoy_a_bemio_output.vtp')
# data = mesh.read(file_name='./data/buoy_a_translate_0p72.vtp')
buoy_a = mesh.read(file_name='./data/buoy_a.GDF', translate=np.array([0.,0.,-.72]))
buoy_a_cut = mesh.cut_mesh(buoy_a, plane_loc=-1e-5)
buoy_a_cut.write(mesh_format='NEMOH')

buoy_b = mesh.read(file_name='./data/buoy_b.GDF', translate=np.array([0.,0.,-21.29]))
buoy_b_cut = mesh.cut_mesh(buoy_b, plane_loc=-1e-5)
buoy_b_cut.write(mesh_format='NEMOH')


# data = mesh.read(file_name='./data/ellipsoid.GDF')
# data.calculate_center_of_gravity()

# Write meshes - uncomment these lines to wrtie the mesh to a different format
# data.write(mesh_format='VTK')
# data.write(mesh_format='NEMOH')
# data.write(mesh_format='WAMIT')

# View mesh - uncomment this line to view the mesh
# data.view()
