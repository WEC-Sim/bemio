#! /usr/bin/python

from bemio.mesh_utilities.mesh import read

flap1 = read(file_name='flap1.gdf')
flap1.view(save_png=False,camera_pos=[70,0,0])
flap1.write(mesh_format='VTP')

flap2 = read(file_name='flap2.gdf')
flap2.view(save_png=False,camera_pos=[70,0,0])
flap2.write(mesh_format='VTP')

base = read(file_name='base.gdf')
base.view(save_png=False,camera_pos=[70,0,0])
base.write(mesh_format='VTP')
