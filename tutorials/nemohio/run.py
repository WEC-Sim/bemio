# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 11:54:47 2014

@author: mlawson
"""
import nemohio as nio
import meshio as mio

mesh = mio.readVtp('/Users/mlawson/Applications/nemoh/matlabRoutines/nonsymmetrical-osx/NonSymmetrical.vtp')
nemoh = nio.Nemoh(simDir='/Users/mlawson/Applications/nemoh/matlabRoutines/nonsymmetrical-py')
nemoh.mesh = mesh
#py.waveFreq = [41,      0.1,     2.0]
#py.waterDepth = 0
#py.runNemohPreProc()
#py.runNemoh()
#py.runNemohPostProc()
nemoh.results.readCMCA()
nemoh.results.plotAddedMassAndDamping()
