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

Author: Michael Lawson
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
    def __init__(self, out_file, density=1000., gravity=9.81):

        self.files = hd.generate_file_names(out_file)
        
        self.rho = density
        self.g = gravity

        self.data = {}
        self._read()


    def _read(self):
        '''
        Function to read WAMIT output file into the class
        Inputs: None
        Outputs: None
        '''
        with open(self.files['out'],'r') as fid:

            raw = fid.readlines()
   
        code = 'WAMIT'
        num_bodies = 0 # Total number of bodies
        bodCount = 0 # Counter for bodies
        freqCount = 0
        T = []
        cg = {}
        cb = {}
        name = {}    
        disp_vol = {}
        k = {}

        
        for i, line in enumerate(raw):
            
            if 'Water depth:' in line:
                water_depth = raw[i].split()[2]
                try:
                    water_depth = np.float(water_depth)
                except:
                    pass 

            # If there is one body in the WAMIT run
            if "Input from Geometric Data File:" in line:

                num_bodies = 1
                name[0] = raw[i].split()[-1]


            
            # If there are two bodies in the WAMIT run
            if "Input from Geometric Data Files:" in line:

                for j in xrange(20): # look for bodies within the next 20 lines

                    if "N=" in raw[i+j]:

                        num_bodies += 1
                        name[num_bodies-1] = raw[i+j].split()[-1]



            # Read the body positions
            if "Total panels:" in line:

                for j in xrange(15): # look for position within the next 15 lines - will only work for wamit files of about 5 bodies

                    if 'XBODY =' in raw[i+j]:
                        '''
                        Note that this is the XBOD YBOD ZBOD defined in the wamit .out file, not the cg as defined in the wamit file
                        '''

                        temp = raw[i+j].split()
                        cg[bodCount] = np.array([temp[2],temp[5],temp[8]]).astype(float)
                        
                    if 'Volumes (VOLX,VOLY,VOLZ):' in raw[i+j]:

                        temp = raw[i+j].split()
                        disp_vol[bodCount] = float(temp[-1])
                        
                    if 'Center of Buoyancy (Xb,Yb,Zb):' in raw[i+j]:

                        temp = raw[i+j].split()
                        cb[bodCount] = np.array([temp[-3],temp[-2],temp[-1]]).astype(float)
                        
                    if 'C(3,3),C(3,4),C(3,5):' in raw[i+j]:

                        temp = np.zeros([6,6])
                        temp2 = raw[i+j].split()
                        temp[2,2] = np.float(temp2[1])
                        temp[2,3] = np.float(temp2[2])
                        temp[2,4] = np.float(temp2[3])
                        
                        temp2 = raw[i+j+1].split()
                        temp[3,3] = np.float(temp2[1])
                        temp[3,4] = np.float(temp2[2])
                        temp[3,5] = np.float(temp2[3])
                        
                        temp2 = raw[i+j+2].split()
                        temp[4,4] = np.float(temp2[1])
                        temp[4,5] = np.float(temp2[2])
                        
                        k[bodCount] = temp
                        
                        
                bodCount += 1      


                
            # Inf freq added mass
            if "Wave period = zero" in line:
                
                amInf  = raw[i+7:i+7+(6*num_bodies)**2]
                amInf = np.array([amInf[temp].split()[2] for temp in xrange(np.size(amInf))]).astype(float)
                amInf = amInf.reshape(6*num_bodies,6*num_bodies)


                
            # Zero freq added mass
            if "Wave period = infinite" in line:
                
                amZero = raw[i+7:i+7+(6*num_bodies)**2]
                amZero = np.array([amZero[temp].split()[2] for temp in xrange(np.size(amZero))]).astype(float)
                amZero = amZero.reshape(6*num_bodies,6*num_bodies)

            # Added mass, damping, and excitation
            if "Wave period (sec)" in line:
                
                T.append(raw[i].split()[4])
                
                am = raw[i+7:i+7+(6*num_bodies)**2]
                am = np.array([am[temp].split()[2] for temp in xrange(np.size(am))]).astype(float)
                am = am.reshape(6*num_bodies,6*num_bodies,1)
                
                rad = raw[i+7:i+7+(6*num_bodies)**2]
                rad = np.array([rad[temp].split()[3] for temp in xrange(np.size(rad))]).astype(float)
                rad = rad.reshape(6*num_bodies,6*num_bodies,1)
                
                ex = raw[i+17+(6*num_bodies)**2:i+17+(6*num_bodies)**2+6*num_bodies]
                ex = np.array([ex[temp].split()[1] for temp in xrange(np.size(ex))]).astype(float)
                ex = ex.reshape(1,6*num_bodies)  
                
                phase = raw[i+17+(6*num_bodies)**2:i+17+(6*num_bodies)**2+6*num_bodies]
                phase = np.array([phase[temp].split()[2] for temp in xrange(np.size(phase))]).astype(float)
                phase = phase.reshape(1,6*num_bodies)  

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
        print 'Dimensionalized WAMIT Hydrodynamic coefficients with g = ' + str(self.g) + ' and rho = ' + str(self.rho)    
        for i in xrange(num_bodies):       
            self.data[i] = hd.HydrodynamicData() 
            self.data[i].name = name[i][0:-4]
            self.data[i].g = self.g
            self.data[i].water_depth = water_depth
            self.data[i].rho = self.rho            
            self.data[i].num_bodies = num_bodies
            self.data[i].body_num = i
            self.data[i].cg = cg[i] 
            self.data[i].cb = cb[i]
            self.data[i].k = k[i]*self.rho*self.g
            self.data[i].disp_vol = disp_vol[i]
            
            self.data[i].am.inf = amInf[6*i:6+6*i,:]
            self.data[i].am.inf = self.data[i].am.inf*self.rho

            self.data[i].am.zero = amZero[6*i:6+6*i,:]
            self.data[i].am.zero = self.data[i].am.zero*self.rho
            
            self.data[i].T = T
            self.data[i].w = 2.0*np.pi/self.data[i].T

            self.data[i].am.all = amAll[6*i:6+6*i,:,:]
            self.data[i].am.all = self.data[i].am.all*self.rho
            
            self.data[i].rd.all = radAll[6*i:6+6*i,:,:]
            for j in xrange(np.shape(self.data[i].rd.all)[2]):
                self.data[i].rd.all[:,:,j] = self.data[i].rd.all[:,:,j]*self.rho*self.data[i].w[j]
                
            self.data[i].ex.mag = exAll[:,6*i:6+6*i]*self.rho*self.g
            self.data[i].ex.phase = np.deg2rad(phaseAll[:,6*i:6+6*i])
            self.data[i].ex.re = self.data[i].ex.mag*np.cos(self.data[i].ex.phase)
            self.data[i].ex.im = self.data[i].ex.mag*np.sin(self.data[i].ex.phase)

            self.data[i].bem_raw_data = raw
            self.data[i].bem_code = code


def clean(directory='.'):

    if _platform == 'darwin' or _platform == 'linux':

        _sys('rm -rf *.out *.1 *.3 wamitlog.txt errorf.log *.p2f *.hst *.p *.h5')

    else:
        print 'The clean function is only supported for osx and linux'


