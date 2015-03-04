from bemio.io import mesh as mio

# Read mesh
mesh = mio.readVtp('./data/non-symmetrical.vtp')

# Write meshes
mesh.writeVtp('./data/non-symmetrical_example_output.vtp')
mesh.writeNemoh('./data/non-symmetrical_example_output.dat')
mesh.writeGdf('./data/non-symmetrical_example_output.gdf')