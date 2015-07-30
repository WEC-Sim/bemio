#! /usr/bin/python

from bemio.mesh_utilities.mesh import read

mesh = read(file_name='coer_comp.gdf')
mesh.view(save_png=False,camera_pos=[2,2,2])
mesh.write(mesh_format='VTP')
