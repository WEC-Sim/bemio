#! /usr/bin/python

from bemio.mesh_utilities.mesh import read

buoy_a = read(file_name='buoy_a.GDF')
buoy_a.view(save_png=False,camera_pos=[50,50,0])
buoy_a.write(mesh_format='VTP')

buoy_b = read(file_name='buoy_b.GDF')
buoy_b.view(save_png=False,camera_pos=[50,50,0])
buoy_b.write(mesh_format='VTP')
