#! /usr/bin/python

from bemio.mesh_utilities.mesh import read

mesh = read(file_name='ellipsoid.gdf')
mesh.view(save_png=False,camera_pos=[20,0,0])
mesh.write(mesh_format='VTP')
