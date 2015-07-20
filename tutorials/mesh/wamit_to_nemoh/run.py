# Import the bemio.mesh_utilities module
from bemio.mesh_utilities import mesh
import numpy as np

# Read WAMIT mesh
buoy = mesh.read(file_name='buoy.gdf')

# Save to a NEMOH mesh
buoy.write(mesh_format='NEMOH')
