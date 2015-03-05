"""
Copyright 2014 the National Renewable Energy Laboratory

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
from bemio.data import bem as hd
from os import system as _sys
from sys import platform as _platform

class WamitOutput(object):
    '''
    Class to read and interact with WAMIT simulation data
    
    Inputs:
        direcotry: location of the wamit output data
        simName: base name of the wamit simulation files
    Outputs:
        None
    '''
    def __init__(self,outFile):

        self.files = hd.generateFileNames(outFile)
        
        self.density = 1000.

        self.data = {}
        self._readOutFile()


    
    def _readOutFile(self):
        '''
        Function to read WAMIT output file into the class
        Inputs: None
        Outputs: None
        '''
        with open(self.files['out'],'r') as fid:

            wamitOut = fid.readlines()
   
        nBodies = 0 # Total number of bodies
        bodCount = 0 # Counter for bodies
        freqCount = 0
        T = []
        cg = {}
        cb = {}
        name = {}    
        volDisp = {}
        k = {}
        pos = {}
        
        for i, line in enumerate(wamitOut):

            # Read gravity and density
            if 'Gravity:' in line:
                gravity = wamitOut[i].split()[1]
                gravity = np.float(gravity)
            
            if 'Water depth:' in line:
                waterDepth = wamitOut[i].split()[2]
                try:
                    waterDepth = np.float(waterDepth)
                except:
                    pass 

            # If there is one body in the WAMIT run
            if "Input from Geometric Data File:" in line:

                nBodies = 1
                name[0] = wamitOut[i].split()[-1]


            
            # If there are two bodies in the WAMIT run
            if "Input from Geometric Data Files:" in line:

                for j in xrange(20): # look for bodies within the next 20 lines

                    if "N=" in wamitOut[i+j]:

                        nBodies += 1
                        name[nBodies-1] = wamitOut[i+j].split()[-1]



            # Read the body positions
            if "Total panels:" in line:

                for j in xrange(20): # look for position within the next 20 lines - will only work for wamit files of about 5 bodies

                    if 'XBODY =' in wamitOut[i+j]:

                        temp = wamitOut[i+j].split()
                        pos[bodCount] = np.array([temp[2],temp[5],temp[8]]).astype(float)
                        
                    if 'Volumes (VOLX,VOLY,VOLZ):' in wamitOut[i+j]:

                        temp = wamitOut[i+j].split()
                        volDisp[bodCount] = np.array([temp[-3],temp[-2],temp[-1]]).astype(float)
                        
                    if 'Center of Buoyancy (Xb,Yb,Zb):' in wamitOut[i+j]:

                        temp = wamitOut[i+j].split()
                        cb[bodCount] = np.array([temp[-3],temp[-2],temp[-1]]).astype(float)
                        
                    if 'C(3,3),C(3,4),C(3,5):' in wamitOut[i+j]:

                        temp = np.zeros([6,6])
                        temp2 = wamitOut[i+j].split()
                        temp[2,2] = np.float(temp2[1])
                        temp[2,3] = np.float(temp2[2])
                        temp[2,4] = np.float(temp2[3])
                        
                        temp2 = wamitOut[i+j+1].split()
                        temp[3,3] = np.float(temp2[1])
                        temp[3,4] = np.float(temp2[2])
                        temp[3,5] = np.float(temp2[3])
                        
                        temp2 = wamitOut[i+j+2].split()
                        temp[4,4] = np.float(temp2[1])
                        temp[4,5] = np.float(temp2[2])
                        
                        k[bodCount] = temp
                        
                    if 'Center of Gravity  (Xg,Yg,Zg):' in wamitOut[i+j]:
                            
                        temp = wamitOut[i+j].split()
                        cg[bodCount] = np.array([temp[-3],temp[-2],temp[-1]]).astype(float)                                                
                        del temp
                        
                        
                bodCount += 1      


                
            # Inf freq added mass
            if "Wave period = zero" in line:
                
                amInf  = wamitOut[i+7:i+7+(6*nBodies)**2]
                amInf = np.array([amInf[temp].split()[2] for temp in xrange(np.size(amInf))]).astype(float)
                amInf = amInf.reshape(6*nBodies,6*nBodies)


                
            # Zero freq added mass
            if "Wave period = infinite" in line:
                
                amZero = wamitOut[i+7:i+7+(6*nBodies)**2]
                amZero = np.array([amZero[temp].split()[2] for temp in xrange(np.size(amZero))]).astype(float)
                amZero = amZero.reshape(6*nBodies,6*nBodies)

            # Added mass, damping, and excitation
            if "Wave period (sec)" in line:
                
                T.append(wamitOut[i].split()[4])
                
                am = wamitOut[i+7:i+7+(6*nBodies)**2]
                am = np.array([am[temp].split()[2] for temp in xrange(np.size(am))]).astype(float)
                am = am.reshape(6*nBodies,6*nBodies,1)
                
                rad = wamitOut[i+7:i+7+(6*nBodies)**2]
                rad = np.array([rad[temp].split()[3] for temp in xrange(np.size(rad))]).astype(float)
                rad = rad.reshape(6*nBodies,6*nBodies,1)
                
                ex = wamitOut[i+17+(6*nBodies)**2:i+17+(6*nBodies)**2+6*nBodies]
                ex = np.array([ex[temp].split()[1] for temp in xrange(np.size(ex))]).astype(float)
                ex = ex.reshape(1,6*nBodies)  
                
                phase = wamitOut[i+17+(6*nBodies)**2:i+17+(6*nBodies)**2+6*nBodies]
                phase = np.array([phase[temp].split()[2] for temp in xrange(np.size(phase))]).astype(float)
                phase = phase.reshape(1,6*nBodies)  

                if freqCount is 0:

                    amAll = am
                    radAll = rad
                    exAll = ex
                    phaseAll = phase
                    
                    freqCount = 1

                else:

                    amAll = np.append(amAll,am,axis=2)
                    radAll = np.append(radAll,rad,axis=2)
                    exAll = np.append(exAll,ex,axis=0)
                    phaseAll = np.append(phaseAll,phase,axis=0)
                                        
        T = np.array(T).astype(float)


        # Load data into the hydrodata strucuture          
        for i in xrange(nBodies):       
            self.data[i] = hd.HydrodynamicData() 
            self.data[i].name = name[i][0:-4]
            self.data[i].g = gravity
            self.data[i].waterDepth = waterDepth
            self.data[i].rho = self.density            
            self.data[i].nBodies = nBodies
            self.data[i].bodyN = i
            self.data[i].cg = cg[i]
            self.data[i].cb = cb[i]
            self.data[i].k = k[i]*self.data[i].rho*self.data[i].g
            self.data[i].pos = pos[i]
            self.data[i].volDisp = volDisp[i]
            
            self.data[i].am.infFreq = amInf[6*i:6+6*i,:]
            self.data[i].am.infFreq = self.data[i].am.infFreq*self.density

            self.data[i].am.zeroFreq = amZero[6*i:6+6*i,:]
            self.data[i].am.zeroFreq = self.data[i].am.zeroFreq*self.density
            
            self.data[i].T = T
            self.data[i].w = 2.0*np.pi/self.data[i].T

            self.data[i].am.all = amAll[6*i:6+6*i,:,:]
            self.data[i].am.all = self.data[i].am.all*self.density
            
            self.data[i].rd.all = radAll[6*i:6+6*i,:,:]
            for j in xrange(np.shape(self.data[i].rd.all)[2]):
                self.data[i].rd.all[:,:,j] = self.data[i].rd.all[:,:,j]*self.density*self.data[i].w[j]
                
            self.data[i].ex.mag = exAll[:,6*i:6+6*i]*self.data[i].rho*self.data[i].g
            self.data[i].ex.phase = np.deg2rad(phaseAll[:,6*i:6+6*i])
            self.data[i].ex.re = self.data[i].ex.mag*np.cos(self.data[i].ex.phase)
            self.data[i].ex.im = self.data[i].ex.mag*np.sin(self.data[i].ex.phase)


def clean(directory='.'):

    if _platform == 'darwin' or _platform == 'linux':

        _sys('rm -rf *.out *.1 *.3 wamitlog.txt errorf.log *.p2f *.hst *.p *.h5')

    else:
        print 'The clean function is only supported for osx and linux'


