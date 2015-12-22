"""
Copyright 2014 the National Renewable Energy Laboratory and Sandia Corporation

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np

import os

from shutil import rmtree as rmd
    
class Nemoh(object):
    '''
    Class that allows simplificed interaction with the Nemoh code
    
    Inputs:
        simDir: The base directory for time Nemoh simulation
    '''
    def __init__(self,simDir):        
        self.dir = simDir
        self.meshDir = os.path.join(self.dir,'mesh')
        self.resultsDir = os.path.join(self.dir,'results')
        if os.path.exists(self.meshDir) is False:
            os.mkdir(self.meshDir)
        if os.path.exists(self.resultsDir) is False:
            os.mkdir(self.resultsDir)
            
        self.gravity = 9.81
        self.density = 1000
        self.name = '.'
        
        self._wavePeriod = None
        self._mesh = None
        
        # These should not be hardcoded
        self.nemohPreProc    = '/Users/mlawson/bin/nemohPreProc'
        self.nemoh              = '/Users/mlawson/bin/nemoh'
        self.nemohPostProc      = '/Users/mlawson/bin/nemohPostProc'
        
        self.files = {}
        self.files['ID.dat']                    = os.path.join(self.dir,'ID.dat')
        self.files['input.txt']                 = os.path.join(self.dir,'input.txt')
        self.files['Nemoh.cal']                 = os.path.join(self.dir,'Nemoh.cal')
        self.files['nemohPreProc.log']          = os.path.join(self.dir,'nemohPreProc.log')
        self.files['nemohPostProc.log']         = os.path.join(self.dir,'nemohPostProc.log')
        self.files['nemoh.log']                 = os.path.join(self.dir,'nemoh.log')
        self.files['nemohMesh.log']             = os.path.join(self.dir,'nemohMesh.log')
        self.files['nemohMesh.dat-forNemohCal']     =   None
        self.files['nemohMesh.dat']     =   None
        

        self.results    = NemohResults(self.dir)
        
        self.writeId()
        
    @property 
    def wavePeriod(self):
        return self._wavePeriod
    @wavePeriod.setter
    def wavePeriod(self,value):
        self._wavePeriod = value
        self.waveFreq = [value[0],2*np.pi/value[2],2*np.pi/value[1]]
        
    @property 
    def mesh(self):
        return self._mesh
    @mesh.setter
    def mesh(self,meshObj):
        self._mesh = meshObj
        self.files['nemohMesh.dat']             = os.path.join(self.meshDir,os.path.split(self._mesh.meshFileName)[1][:-3]+'dat')
        self.files['nemohMesh.dat-forNemohCal'] = os.path.join('mesh',os.path.split(self._mesh.meshFileName)[1][:-3]+'dat')
        self._mesh.writeNemohMesh(self.files['nemohMesh.dat'])
    
    def clean(self):
        '''
        Function to remove all file created by Nemoh from the case directory
        Inputs:
            None
        Otuputs:
            None
        '''
        rmd(self.dir)
        os.system('rm ' + self.baseDir + os.path.sep + '*.log')
        try:
            os.remove(self.files['ID.dat'])
        except:
            pass
        try:
            os.remove(self.mesh.files['Mesh.cal'])
        except:
            pass
                    
    def runNemohPreProc(self):
        '''
        Function to run the Nemoh Pre Processor
        Inputs:
            None
        Outputs:
            None
        '''
        os.chdir(self.dir)
        self.writeNemohCal()
        if os.sys.platform == 'darwin':
            os.system(self.nemohPreProc + ' | tee ' + self.files['nemohPreProc.log'])
        else:
            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
            
    def runNemoh(self):
        '''
        Function to run the Nemoh
        Inputs:
            None
        Outputs:
            None
        '''
        self.writeInput()
        self.writeNemohCal()
        os.chdir(self.dir)
        if os.sys.platform == 'darwin':
            os.system(self.nemoh + ' | tee ' + self.files['nemoh.log'])
        else:
            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
            
    def runNemohPostProc(self):
        '''
        Function to run the Nemoh Post Processor
        Inputs:
            None
        Outputs:
            None
        '''
        self.writeNemohCal()
        os.chdir(self.dir)
        if os.sys.platform == 'darwin':
            os.system(self.nemohPostProc + ' | tee ' + self.files['nemohPostProc.log'])
        else:
            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
            
    def writeId(self):
        '''
        Function to write the ID.dat file for Neoh simulations
        Inputs:
            None
        Outputs:
            None
        '''
        with open(self.files['ID.dat'],'w') as fid:
            fid.write(str(len(self.name)))
            fid.write('\n')
            fid.write(self.name)  

    def writeInput(self):
        '''
        Function to write the inout.txt file for Neoh simulations
        Inputs:
            None
        Outputs:
            None
        '''
        with open(self.files['input.txt'],'w') as fid:
            fid.write('\n0')
            
    def writeNemohCal(self):
        nemohCalFile = []
        nemohCalFile.append('--- Environment ------------------------------------------------------------------------------------------------------------------ \n')
        nemohCalFile.append(str(self.density) + ' ! RHO ! KG/M**3 ! Fluid specific volume \n')
        nemohCalFile.append(str(self.gravity) + '                            ! G                     ! M/S**2        ! Gravity \n')
        nemohCalFile.append(str(self.waterDepth) + '                ! DEPTH                       ! M             ! Water depth\n')
        nemohCalFile.append('0.      0.              ! XEFF YEFF             ! M             ! Wave measurement point\n')
        nemohCalFile.append('--- Description of floating bodies -----------------------------------------------------------------------------------------------\n')
        nemohCalFile.append('1                               ! Number of bodies\n')
        nemohCalFile.append('--- Body 1 -----------------------------------------------------------------------------------------------------------------------\n')
        nemohCalFile.append("'" + self.files['nemohMesh.dat-forNemohCal'] + "'" + '             ! Name of mesh file\n')
        nemohCalFile.append(str(self._mesh.nPoints) + ' ' + str(self._mesh.nFaces) + '                  ! Number of points and number of panels         \n')
        nemohCalFile.append('6                               ! Number of degrees of freedom\n')
        nemohCalFile.append('1 1. 0. 0. 0. 0. 0.             ! Surge\n')
        nemohCalFile.append('1 0. 1. 0. 0. 0. 0.             ! Sway\n')
        nemohCalFile.append('1 0. 0. 1. 0. 0. 0.             ! Heave\n')
        nemohCalFile.append('2 1. 0. 0. 0. 0. 0.000000               ! Roll about a point\n')
        nemohCalFile.append('2 0. 1. 0. 0. 0. 0.000000               ! Pitch about a point\n')
        nemohCalFile.append('2 0. 0. 1. 0. 0. 0.000000               ! Yaw about a point\n')
        nemohCalFile.append('6                               ! Number of resulting generalised forces\n')
        nemohCalFile.append('1 1. 0. 0. 0. 0. 0.             ! Force in x direction\n')
        nemohCalFile.append('1 0. 1. 0. 0. 0. 0.             ! Force in y direction\n')
        nemohCalFile.append('1 0. 0. 1. 0. 0. 0.             ! Force in z direction\n')
        nemohCalFile.append('2 1. 0. 0. 0. 0. 0.000000               ! Moment force in x direction about a point\n')
        nemohCalFile.append('2 0. 1. 0. 0. 0. 0.000000               ! Moment force in y direction about a point\n')
        nemohCalFile.append('2 0. 0. 1. 0. 0. 0.000000               ! Moment force in z direction about a point\n')
        nemohCalFile.append('0                               ! Number of lines of additional information \n')
        nemohCalFile.append('--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n')
        nemohCalFile.append(str(self.waveFreq).replace(',','').replace('[','').replace(']','') + '           ! Number of wave frequencies, Min, and Max (rad/s)\n')
        nemohCalFile.append('1       0.000000 0.000000               ! Number of wave directions, Min and Max (degrees)\n')
        nemohCalFile.append('--- Post processing ---------------------------------------------------------------------------------------------------------------\n')
        nemohCalFile.append('0       0.1     10.             ! IRF                           ! IRF calculation (0 for no calculation), time step and duration\n')
        nemohCalFile.append('0                               ! Show pressure\n')
        nemohCalFile.append('0       0.      180.            ! Kochin function               ! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n')
        nemohCalFile.append('0       50      400.    400.    ! Free surface elevation        ! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction\n')
        nemohCalFile.append('-----\n')

        with open(self.files['Nemoh.cal'],'w') as fid:
            fid.writelines(nemohCalFile)

        
#    def writeMeshCal(self):
#        with open(self.files['Mesh.cal'],'w') as fid:
#            fid.write(self.files['Mesh-noPath'])
#            fid.write('\n')
#            fid.write('0')
#            fid.write('\n')
#            fid.write('0. 0.') # not sure what this is
#            fid.write('\n')
#            fid.write(str(self._cg).replace('[','').replace(']','').replace(',',''))
#            fid.write('\n')
#            fid.write(str(self.targetNumPanels))
#            fid.write('\n')
#            fid.write('2')
#            fid.write('\n')
#            fid.write('0')
#            fid.write('\n')
#            fid.write('1')
        
#    def runNemohMesh(self):
#        os.chdir(self.dir)
#        self.writeMeshCal()
#        self.writeNemohMeshInput()
#        if os.sys.platform == 'darwin':
#            os.system(self.nemohMesh + ' | tee ' + self.files['nemohMesh.log'])
#        else:
#            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
#        self.readNemohMesh()