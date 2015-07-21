#! /usr/bin/python

from bemio.mesh_utilities.mesh import read

flap = read(file_name='flap.gdf')
flap.view(save_png=False,camera_pos=[70,0,0])
flap.write(mesh_format='VTP')

base = read(file_name='base.gdf')
base.view(save_png=False,camera_pos=[70,0,0])
base.write(mesh_format='VTP')
