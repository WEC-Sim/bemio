# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 14:36:03 2014

@author: mlawson
"""

"""
Created on Tue Nov  4 15:28:52 2014

@author: mlawson

This module reads data from Nemoh simulations and convert from GDF to Nemoh mesh format
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from shutil import rmtree as rmd
    
class Nemoh(object):
    def __init__(self,directory,name):        
        self.baseDir = directory
        self.name = name
        self.dir = os.path.join(directory, name)
        if os.path.exists(self.baseDir) is False:
            os.mkdir(self.dir)
        del directory, name
        if os.path.exists(self.dir) is False:
            os.mkdir(self.dir)
        
        self.files = {}
        self.files['ID.dat'] = os.path.join(self.baseDir,'ID.dat')

        self.mesh       = NemohMesh(self)
        self.sim        = NemohSimulation(self)   
        self.results    = NemohResults(self)
        
        self.writeId()
        
    def writeId(self):        
        with open(self.files['ID.dat'],'w') as fid:
            fid.write(str(len(self.name)))
            fid.write('\n')
            fid.write(self.name)  
    
    def clean(self):
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

class NemohMesh(object):
    def __init__(self,sim):
        self.dir = os.path.join(sim.dir,'mesh')
        self.baseDir = os.path.join(sim.dir, '..')
        self.name = sim.name
        if os.path.exists(self.dir) is False:
            os.mkdir(self.dir)
        
        self.files = {}
        self.files['NemohMesh.tec']             = os.path.join(self.dir,'Mesh.tec')
        self.files['L12.dat']                   = os.path.join(self.dir,'L12.dat')
        self.files['L10.dat']                   = os.path.join(self.dir,'L10.dat')
        self.files['Kochin.dat']                = os.path.join(self.dir,'Kochin.dat')
        self.files['KH.dat']                    = os.path.join(self.dir,'KH.dat')
        self.files['Integration.dat']           = os.path.join(self.dir,'Integration.dat')
        self.files['Inertia_hull.dat']          = os.path.join(self.dir,'Inertia_hull.dat')
        self.files['Hydrostatics.dat']          = os.path.join(self.dir,'Hydrostatics.dat')
        self.files['CG_hull.dat']               = os.path.join(self.dir,'CG_hull.dat')
        self.files['Freesurface.dat']           = os.path.join(self.dir,'Freesurface.dat')
        self.files['Description_Wetted.tec']    = os.path.join(self.dir,'Description_Wetted.tec')
        self.files['Description_Full.tec']      = os.path.join(self.dir,'Description_Full.tec')
        self.files['Mesh.cal']                  = os.path.join(self.baseDir,'Mesh.cal')

        self.faces = []
        self.cords = []
        
#        self.nemohFaces = []
#        self.nemohCords = []
        
#        self._gdfFile = None
#        self._stlFile = None
#        self._vtpFile = None
        self._cg = None
        self._symmetry = None
        self._targedNumPanels = None
        self._numBodies = None
        
        # These should not be hardcoded
        self.nemohMesh          = '/Users/mlawson/bin/nemohMesh'
        
    def setMeshFileNames(self):
#        self.files['nemohMeshInput'] = self.dir + os.path.sep + self.meshName + 'nemohMesh'
        self.files['nemohMesh'] =               os.path.join(self.dir,self.meshName+'.dat')
        self.files['nemohMesh.log'] =           os.path.join(self.baseDir + os.path.sep + 'nemohMesh.log')
        self.files['nemohMesh-forNemohCal'] =   os.path.join(self.name,'mesh',self.meshName,'nemohMesh.dat')       
        temp, self.files['Mesh-noPath'] =       os.path.split(self.files['nemohMesh'][:-4])
        
#    @property
#    def stlFile(self):
#        return self._stlFile
#    @stlFile.setter
#    def stlFile(self,fName):
#        self._stlFile = fName
#        self.files['stl'] = self._stlFile
#        temp, self.meshName = os.path.split(self._stlFile[:-3])
#        self.setMeshFileNames()
#        self.readSTL()
#        self.writeNemohMeshInputVtp()
#        
#    @property
#    def gdfFile(self):
#        return self._gdfFile
#    @gdfFile.setter
#    def gdfFile(self,fName):
#        self._gdfFile = fName
#        self.files['GDF'] = self._gdfFile
#        temp, self.meshName = os.path.split(self._gdfFile[:-3])
#        self.setMeshFileNames()
#        self.readGDF()
#        self.writeNemohMeshInputVtp()
#        
#    @property
#    def vtpFile(self):
#        return self._vtpFile
#    @vtpFile.setter
#    def vtpFile(self,fName):
#        self._vtpFile = fName
#        self.files['vtp'] = self._vtpFile
#        temp, self.meshName = os.path.split(self._vtpFile[:-3])
#        self.setMeshFileNames()
#        self.readVTP()
#        self.writeNemohMeshInputVtp()
        
    @property
    def cg(self):
        return self._cg        
    @cg.setter
    def cg(self,cg):
        if np.shape(cg) == (3,) is False:
            raise Exception('cg input error')
        self._cg = cg
        
    @property
    def symmetry(self):
        return self._symmetry        
    @symmetry.setter
    def symmetry(self,symmetry):
        if np.shape(symmetry) == (3,) is False:
            raise Exception('symmetry input error')
        self._symmetry = symmetry
        
    @property 
    def targedNumPanels(self):
        return self
    @targedNumPanels.setter
    def targedNumPanels(self,numPanels):
        self._targedNumPanels = numPanels
        if type(self._targedNumPanels) is not int:
            raise Exception('targedNumPanels must be an integer')
            
    @property 
    def numBodies(self):
        return self
    @numBodies.setter
    def numBodies(self,numBodies):
        self._numBodies = numBodies
        if type(self._numBodies) is not int:
            raise Exception('numBodies must be an integer')    
        if self._numBodies != 1:
            raise Exception('numBodies must be 1')  
        
    def removeSurfacePanels(self):
        tempFaces = []
        count = 0
        
        for i in xrange(self.numFaces):
            deleteFace = 0
            p0 = self.faces[i][0]
            p1 = self.faces[i][1]
            p2 = self.faces[i][2]
            p3 = self.faces[i][3]
            z0 = float(self.cords[int(p0)][2])
            z1 = float(self.cords[int(p1)][2])
            z2 = float(self.cords[int(p2)][2])
            z3 = float(self.cords[int(p3)][2])
            
            if z0 == 0.:
                deleteFace += 1
            if z1 == 0.:
                deleteFace += 1
            if z2 == 0.:
                deleteFace += 1
            if z3 == 0.:
                deleteFace += 1
            if deleteFace != 4:
                tempFaces.append(self.faces[i])
                count  += 1
        print 'removed ' + str(count) + ' surface faces'
        self.faces = []
        self.faces = tempFaces
        self.numFaces = np.shape(self.faces)[0]
        os.remove(self.files['MeshInputVtp'])
        self.writeNemohMeshInputVtp()
        
    def writeNemoMeshOutputVtp(self):
        writeVtp(self.nemohCords,self.nemohFaces,self.files['MeshOutputVtp'])
        
    def writeNemohMeshInputVtp(self):
        writeVtp(self.cords,self.faces,self.files['MeshInputVtp'])
                
#    def writeNemohMeshInput(self):   
#        with open(self.files['nemohMeshInput'],'w') as fid:
#            fid.write(str(self.numCords))
#            fid.write('\n')
#            fid.write(str(self.numFaces))
#            fid.write('\n')
#            for i in range(np.shape(self.cords)[0]):
#                fid.writelines(str(self.cords[i]).replace('[','').replace(']',''))
#                fid.write('\n')
#            for i in range(self.numFaces):
#                fid.write(str(self.faces[i]+1).replace('[','').replace(']',''))
#                fid.write('\n')

    def writeMeshCal(self):
        with open(self.files['Mesh.cal'],'w') as fid:
            fid.write(self.files['Mesh-noPath'])
            fid.write('\n')
            fid.write('0')
            fid.write('\n')
            fid.write('0. 0.') # not sure what this is
            fid.write('\n')
            fid.write(str(self._cg).replace('[','').replace(']','').replace(',',''))
            fid.write('\n')
            fid.write(str(self.targetNumPanels))
            fid.write('\n')
            fid.write('2')
            fid.write('\n')
            fid.write('0')
            fid.write('\n')
            fid.write('1')
            
    def runNemohMesh(self):
        os.chdir(self.baseDir)
        self.writeMeshCal()
        self.writeNemohMeshInput()
        if os.sys.platform == 'darwin':
            os.system(self.nemohMesh + ' | tee ' + self.files['nemohMesh.log'])
        else:
            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
        self.readNemohMesh()
        self.writeNemoMeshOutputVtp()
         
class NemohSimulation(object):

    def __init__(self,sim):
        self.gravity = 9.81
        self.density = 1000.
        self.baseDir = sim.baseDir
        self.dir = sim.dir

        self.files = {}
        self.files['input.txt'] = sim.dir + os.path.sep + 'input.txt'
        self.files['Nemoh.cal'] = sim.dir + os.path.sep + 'Nemoh.cal'
        self.files['nemohPreProc.log'] = sim.baseDir + os.path.sep + 'nemohPreProc.log'
        self.files['nemohPostProc.log'] = sim.baseDir + os.path.sep + 'nemohPostProc.log'
        self.files['nemoh.log'] = sim.baseDir + os.path.sep + 'nemoh.log'
        self.files['Normalvelocities.dat'] = sim.dir + os.path.sep + 'Normalvelocities.dat'
        
        self._wavePeriod = None
        
        # These should not be hardcoded
        self.nemohPreProc    = '/Users/mlawson/bin/nemohPreProc'
        self.nemoh              = '/Users/mlawson/bin/nemoh'
        self.nemohPostProc      = '/Users/mlawson/bin/nemohPostProc'
        
    @property 
    def wavePeriod(self):
        return self
    @wavePeriod.setter
    def wavePeriod(self,value):
        self._wavePeriod = value
        self.waveFreq = [value[0],2*np.pi/value[2],2*np.pi/value[1]]
        
    def writeInput(self):
        with open(self.files['input.txt'],'w') as fid:
            fid.write('\n0')
            
    def writeNemohCal(self,sim):
        nemohCalFile = []
        nemohCalFile.append('--- Environment ------------------------------------------------------------------------------------------------------------------ \n')
        nemohCalFile.append(str(self.density) + ' ! RHO ! KG/M**3 ! Fluid specific volume \n')
        nemohCalFile.append(str(self.gravity) + '                            ! G                     ! M/S**2        ! Gravity \n')
        nemohCalFile.append('10.000000                 ! DEPTH                       ! M             ! Water depth\n')
        nemohCalFile.append('0.      0.              ! XEFF YEFF             ! M             ! Wave measurement point\n')
        nemohCalFile.append('--- Description of floating bodies -----------------------------------------------------------------------------------------------\n')
        nemohCalFile.append('1                               ! Number of bodies\n')
        nemohCalFile.append('--- Body 1 -----------------------------------------------------------------------------------------------------------------------\n')
        nemohCalFile.append(sim.mesh.files['nemohMesh-forNemohCal'] + '             ! Name of mesh file\n')
        nemohCalFile.append(str(sim.mesh.numCords) + ' ' + str(sim.mesh.numFaces) + '                  ! Number of points and number of panels         \n')
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
            
    def runNemohPreProc(self):
        os.chdir(self.baseDir)
        if os.sys.platform == 'darwin':
            os.system(self.nemohPreProc + ' | tee ' + self.files['nemohPreProc.log'])
        else:
            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
            
    def runNemoh(self):
        os.chdir(self.baseDir)
        if os.sys.platform == 'darwin':
            os.system(self.nemoh + ' | tee ' + self.files['nemoh.log'])
        else:
            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
            
    def runNemohPostProc(self):
        os.chdir(self.baseDir)
        if os.sys.platform == 'darwin':
            os.system(self.nemohPostProc + ' | tee ' + self.files['nemohPostProc.log'])
        else:
            raise Exception('Only osx is supported at this time (linux and windows will prob work too with small change to this python module)')
        
class NemohResults(object):
    def __init__(self,sim):
        self.baseDir = sim.baseDir
        self.dir = sim.dir + os.path.sep + 'results'
        if os.path.exists(self.dir) is False:
            os.mkdir(self.dir)
        
        self.files = {}
        self.w = [] # need to clean up so this is consistant with the sim class data
        self.addedMass = {}
        self.radiationDamping = {}
        
        self.files['RadiationCoefficients.tec'] = self.dir + os.path.sep + 'RadiationCoefficients.tec'
        self.files['index.dat'] = self.dir + os.path.sep + 'index.dat'
        self.files['Force.dat'] = self.dir + os.path.sep + 'Force.dat'
        self.files['FKForces.dat'] = self.dir + os.path.sep + 'FKForces.dat'
        self.files['Fe.dat'] = self.dir + os.path.sep + 'Fe.dat'
        self.files['CM.dat'] = self.dir + os.path.sep + 'CM.dat'
        self.files['CA.dat'] = self.dir + os.path.sep + 'CA.dat'
        self.files['Normalvelocities.dat'] = sim.dir + os.path.sep + 'Normalvelocities.dat'
        
    def readCMCA(self):
        '''
        Function to read CM.dat created by Nemoh
        '''        
        # load  data files
        with open(self.files['CM.dat']) as fid :
            linesCM = fid.readlines()        
        with open(self.files['CA.dat']) as fid :
            linesCA = fid.readlines()

        # Read the number of frequencies
        self.numFreqs = int(linesCM[0].split()[-1])

        # Read the Frequencies, the added mass matrix, and the radiation damping matrix at each frequency
        for i in xrange(self.numFreqs):
            self.w.append(float(linesCM[1+i*7].replace('\n','')))
            self.addedMass[self.w[i]] = [temp.replace('\n','') for temp in linesCM[2+i*7:8+i*7]]
            self.addedMass[self.w[i]] = np.array([np.fromstring(temp,sep=' ') for temp in self.addedMass[self.w[i]]])        
            self.radiationDamping[self.w[i]] = [temp.replace('\n','') for temp in linesCA[2+i*7:8+i*7]]
            self.radiationDamping[self.w[i]] = np.array([np.fromstring(temp,sep=' ') for temp in self.radiationDamping[self.w[i]]])
        
        addedMassDiag = []
        radiationDampingDiag = []
        for i in xrange(self.numFreqs):
            
            addedMassDiag.append(np.diag(self.addedMass[self.w[i]]))
            radiationDampingDiag.append(np.diag(self.radiationDamping[self.w[i]]))
            
        self.addedMassDiag = np.array(addedMassDiag)
        self.radiationDampingDiag = np.array(radiationDampingDiag)
             
    def plotAddedMassAndDamping(self):
        f, ax = plt.subplots(2, sharex=True)
        ax[0].plot()
        ax[1].plot()
        
        for i in xrange(6):
            ax[0].set_title('Diagional Compinent of Added Mass Matrix')
            ax[0].plot(self.w,self.addedMassDiag[:,i],'x-',label='Component (' + str(i) + ', ' + str(i) + ')')
            ax[0].set_ylabel('Added Mass (kg)')
            
            ax[1].plot(self.w,self.radiationDampingDiag[:,i],'x-',label='Component (' + str(i) + ', ' + str(i) + ')')
            ax[1].set_title('Diagional Compinent of Radiation Damping Matrix')
            ax[1].set_xlabel('Wave Frequency (rad/s)')
            ax[1].set_ylabel('Radiation Damping')
            ax[1].legend(loc=0)
            
        plt.show()