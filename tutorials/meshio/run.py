import meshio as mio

# Read mesh
mesh = mio.readVtp('./data/non-symmetrical.vtp')

# Write meshes
#mesh.writeVtp('./NonSymmetrical-demoOutput.vtp')
#mesh.writeNemohMesh('./NonSymmetrical-demoOutput.dat')
#mesh.writeGdf('./NonSymmetrical-demoOutput.gdf')