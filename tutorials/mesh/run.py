from bemio.mesh_utilities import mesh 
import matplotlib.pyplot as plt

plt.interactive(True)

# Read mesh
non_symmetrical = mesh.read(file_name='./data/non_symmetrical.gdf')

# Write meshes
non_symmetrical.write_vtp(out_file='./data/non_symmetrical_example_output.vtp')
non_symmetrical.write_nemoh(out_file='./data/non_symmetrical_example_output.dat')
non_symmetrical.write_gdf(out_file='./data/non_symmetrical_example_output.gdf')

# Collapse mesh surface
#non_symmetrical.collapse(value=6.0)

# View mesh
#non_symmetrical.view()