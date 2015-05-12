from bemio.mesh_utilities import mesh 

# Read mesh
non_symmetrical = mesh.read(file_name='./data/non_symmetrical.vtp')

# Write meshes
non_symmetrical.write_vtp(out_file='./data/non_symmetrical_example_output.vtp')
non_symmetrical.write_nemoh(out_file='./data/non_symmetrical_example_output.dat')
non_symmetrical.write_gdf(out_file='./data/non_symmetrical_example_output.gdf')