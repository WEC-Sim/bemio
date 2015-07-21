#! /usr/bin/python

from bemio.mesh_utilities.mesh import read

sphere = read(file_name='sphere.gdf')
sphere.view(save_png=False,camera_pos=[70,0,0])
sphere.write(mesh_format='VTP')
