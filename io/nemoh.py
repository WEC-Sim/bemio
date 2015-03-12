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

This class contains a structure to store hydrodynamic data from WAMTI,
AQWA, Nemoh, or another code that calculates hydrodynamic coefficinets
and excitation forces

Author: Michael Lawson
"""

import os

import numpy as np

from bemio.data import bem as hd

try:

    from astropy.io import ascii
 
except:

    raise Exception('The astropy module must be installed. Try "pip install astropy"')   


class NemohOutput   (object):
    '''
    Class that is used to read Nemoh results.

    Inputs:
    results_dir -- the directory with the Nemoh results files (e.g. CA.dat)
    '''
    def __init__(self, sim_dir='./', cal_file='Nemoh.cal', results_dir = 'Results'):

        self.dir = os.path.abspath(sim_dir)

        self.files = {}        
        self.files['Nemoh']     = os.path.join(self.dir,cal_file)
        # self.files['index']     = os.path.join(self.dir,results_dir,'index.dat')
        # self.files['Force']     = os.path.join(self.dir,results_dir,'Force.dat')
        # self.files['FKForces']  = os.path.join(self.dir,results_dir,'FKForces.dat')
        # self.files['Fe']        = os.path.join(self.dir,results_dir,'Fe.dat')
        # self.files['CM']        = os.path.join(self.dir,results_dir,'CM.dat')
        # self.files['CA']        = os.path.join(self.dir,results_dir,'CA.dat')
        self.files['RadiationCoefficients'] = os.path.join(self.dir,results_dir,'RadiationCoefficients.tec')

        self.data = {}

        self._readRadiationCoefficients()

    def _readRadiationCoefficients(self):
        '''
        Function to read CM.dat and CA.dat files from Nemoh

        Internal function called during at __init__
        '''

        # Read added mass and damping 
        with open(self.files['RadiationCoefficients']) as fid:
            amrd = fid.readlines()

        z_ind = []
        amrdProc = {}  
        for i, line in enumerate(amrd):

            if 'Zone' in line:

                z_ind.append(i)
                
                if np.size(z_ind) > 1:

                    amrdProc[amrd[z_ind[-2]]] = ascii.read(amrd[z_ind[-2]+1:z_ind[-1]])

        amrdProc[amrd[z_ind[-1]]] = ascii.read(amrd[z_ind[-1]+1::]) # read in the last one

        for i, key in enumerate(amrdProc.keys()):

            if i == 0:

                w = np.array(amrdProc[key].field(0))

            am[key] = 










        self.data[0] = hd.HydrodynamicData()
        self.data[0].w = w
        self.data[0].T = 2.*np.pi/w





















        # with open(self.files['CM.dat']) as fid :
        #     cm = fid.readlines()        
        # with open(self.files['CA.dat']) as fid :
        #     ca = fid.readlines()



        # # Read the number of frequencies
        # self.numFreqs = int(cm[0].split()[-1])

        # # Read the Frequencies, the added mass matrix, and the radiation damping matrix at each frequency
        # for i in xrange(self.numFreqs):
        #     self.freq.append(float(cm[1+i*7].replace('\np','')))
        #     self.addedMass[self.freq[i]] = [temp.replace('\n','') for temp in cm[2+i*7:8+i*7]]
        #     self.addedMass[self.freq[i]] = np.array([np.fromstring(temp,sep=' ') for temp in self.addedMass[self.freq[i]]])        
        #     self.radiationDamping[self.freq[i]] = [temp.replace('\n','') for temp in ca[2+i*7:8+i*7]]
        #     self.radiationDamping[self.freq[i]] = np.array([np.fromstring(temp,sep=' ') for temp in self.radiationDamping[self.freq[i]]])
        
        # addedMassDiag = []
        # radiationDampingDiag = []
        # for i in xrange(self.numFreqs):
            
        #     addedMassDiag.append(np.diag(self.addedMass[self.freq[i]]))
        #     radiationDampingDiag.append(np.diag(self.radiationDamping[self.freq[i]]))
            
        # self.addedMassDiag = np.array(addedMassDiag)
        # self.radiationDampingDiag = np.array(radiationDampingDiag)