# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 11:54:47 2014

@author: mlawson
"""
from bemio.io import nemoh as nio


nem = nio.NemohOutput(sim_dir='./data/two_body')
nem.read_hydrostatics(body_num=0,file='./data/two_body/Mesh/Hydrostatics_0.dat')
nem.read_kh(body_num=0, file='./data/two_body/Mesh/KH_0.dat')

