# Import the bemio.mesh_utilities module
from bemio.mesh_utilities import mesh

# Read mesh
non_symmetrical = mesh.read(file_name='./data/non_symmetrical.gdf')

# Write meshes
non_symmetrical.write(mesh_format='VTK')
non_symmetrical.write(mesh_format='NEMOH')
non_symmetrical.write(mesh_format='WAMIT')

# View mesh
non_symmetrical.view()

