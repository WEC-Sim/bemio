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

import numpy as np

import os

import matplotlib.pyplot as plt

from sys import platform as _platform

import pickle

class HydrodynamicCoefficients(object):
    '''
    Data class that contains hydrodynamic coefficient data

    Variables:
    all -- Frequency dependent hydrodynamic coefficients. 
    6 x (6 OR 6*nbodies) x numFrequencies np.array.

    inf -- Infinate freqnency added mass. 6 x (6 OR 6*nbodies)
    np.array()

    zero -- Zero frequency added mass. 6 x (6 OR 6*nbodies) np.array()
    '''
    def __init__(self):
        self.all            = np.array([])
        self.inf            = np.array([])
        self.zero           = np.array([])
    
class HydrodynamicExcitation(object):
    '''
    Data class that contains hydrodynamic excitation coefficinets.

    Variables:
    re -- Real component of hydrodynamic excitation 
    6 x numFrequencies np.array()

    im -- Imaginary component of hydrodynamic excitation 
    6 x numFrequencies np.array()

    mag -- Magnitude of hydrodynamic excitation 
    6 x numFrequencies np.array()

    phase -- Phase angle of hydrodynamic excitation 
    6 x numFrequencies np.array(). Should be in radians

    '''
    def __init__(self):
        self.re             = np.array([])
        self.im             = np.array([])
        self.mag            = np.array([])
        self.phase          = np.array([])
        
class HydrodynamicData(object):
    ''''
    Sturcuture for storing data from WAMIT, AQWA and Nemoh.

    Variables:
    rho -- density
    g -- gravity
    files -- Python dictionary of files associated with the
    simulation
    nBodies -- Number of bodies in simulation
    cg -- Center of gravity
    cb -- Center of buoyancy
    volDisp -- Volume displacement
    T -- Wave periods of simulations (e.g. [1, 2, 3, 4, 5,])
    w -- Wave freqencies of simulations 
    am -- Added mass coefficients of HydrodynamicCoefficients type
    rd -- Radiation damping coefficients of HydrodynamicCoefficients
    type
    wpArea -- Water plane area          
    buoyForce -- Buoyancy force at equilibrium
    k -- Hydrostatic stiffness matrix
    ex -- Excitation coeffs of HydrodynamicExcitation type
    waterDepth -- Water depth
    waveDir -- Wave direction 
    name -- Name of the body in the simulation
    '''

    def __init__(self):
        self.rho            = 1000.
        self.g              = 9.81      
        self.files          = {}
        self.nBodies        = 0                             
        self.bodyN          = None
        self.cg             = {}                            
        self.cb             = {}                            
        self.volDisp        = {}                            
        self.T              = {}                            
        self.w              = {}                            
        self.am             = HydrodynamicCoefficients()    
        self.rd             = HydrodynamicCoefficients()    
        self.wpArea         = {}                            
        self.buoyForce      = {}                            
        self.k              = {}                            
        self.ex             = HydrodynamicExcitation()      
        self.waterDepth     = None                          
        self.waveDir        = 0                             
        self.name           = None                          
        
    def plotAddedMassAndDamping(self,components):
        '''
        Function to plot the added mass and raditation damping coefficinets

        Inputs:
        components -- A list of components to plot. E.g [[0,0],[1,1],[2,2]]
        
        Outputs:
        None -- A plot is displayed. The plt.show() command may need to be used
        depending on your python env settings
        '''                        
        
        f, ax = plt.subplots(2, sharex=True, figsize=(8,10))
        
        # Frame 0 - added mass
        ax[0].plot()
        ax[0].set_title('Hydrodynamic coefficients for body ' + str(self.name))    
        ax[0].set_ylabel('Added mass')
        
        # Frame 1 - radiation damping
        ax[1].plot()
        ax[1].set_xlabel('Wave frequency (rad/s)')
        ax[1].set_ylabel('Radiation damping')
        
        # Plot added mass and damping
        for i,comp in enumerate(components):
            
            x = comp[0]
            y = comp[1]
            w = self.w
            rd = self.rd.all[x,y,:]
            am = self.am.all[x,y,:]

            ax[0].plot(w,am,'x-',label='Component (' + str(x) + ', ' + str(y) + ')')
            ax[1].plot(w,rd,'x-',label='Component (' + str(x) + ', ' + str(y) + ')')
            
        # Show legend on frame 0
        ax[0].legend(loc=0)

    def plotExcitation(self,components):
        '''
        Function to plot wave excitation coefficients
        
        Inputs:
        components -- A list of components to plot. E.g [0,1,2,5]
        
        Outputs:
        None -- A plot is displayed. The plt.show() command may need to be used
        depending on your python env settings
        '''
        
        f, ax = plt.subplots(4, sharex=True,figsize=(8,10))

        # Frame 0 - magnitude
        ax[0].plot()
        ax[0].set_ylabel('Ex force - mag')
        ax[0].set_title('Excitation force for body ' + str(self.name))    

        # Frame 1 - phase
        ax[1].plot()        
        ax[1].set_xlabel('Wave frequency (rad/s)')        
        ax[1].set_ylabel('Ex force - phase')

        # Frame 2 - real
        ax[2].plot()
        ax[2].set_ylabel('Ex force - real')
        
        # Frame 3 - imaginary
        ax[3].plot()
        ax[3].set_ylabel('Ex force - imaginary')
        
        for i,comp in enumerate(components):
            
            m = comp
            w = self.w
            re = self.ex.re[:,m]
            im = self.ex.im[:,m]
            mag = self.ex.mag[:,m]
            phase = self.ex.phase[:,m]

            ax[0].plot(w,mag,'x-',label='Component (' + str(m) + ')')
            ax[1].plot(w,phase,'x-',label='Component (' + str(m) + ')')
            ax[2].plot(w,re,'x-',label='Component (' + str(m) + ')')
            ax[3].plot(w,im,'x-',label='Component (' + str(m) + ')')

            ax[0].legend(loc=0)

def writePickle(data,outFile):
    '''
    Writes hydrodynamic data to a pickle file.
    
    Inputs:
    data -- dictionary that contains HydrodynamicData objects for each body in the simulation
    outFile -- name of the pickle file

    Outputs: None
    '''

    pickle.dump(data,open(outFile,'wb'))
    
    print 'Wrote pickle data to ' + outFile

def writeHdf5(data,outFile):
    '''
    Writes hydrodynamic data to a HDF5 file structure.
    
    Inputs:
    data -- dictionary that contains HydrodynamicData objects for each body in the simulation
    outFile -- name of the hdf5 file

    Outputs: None
    '''
    try:

        import h5py

    except:

        raise Exception('The h5py module must be installed to used the writeHdf5 functionality.')
        
        
        
    with h5py.File(outFile, "w") as f:       

        for key, key in enumerate(data.keys()):

            T = f.create_dataset('body' + str(key) + '/sim/T',data=data[key].T)
            T.attrs['units'] = 's'
            T.attrs['description'] = 'Wave periods'
            
            w = f.create_dataset('body' + str(key) + '/sim/w',data=data[key].w)
            w.attrs['units'] = 'rad/s'                
            w.attrs['description'] = 'Wave frequencies'                
            
            k = f.create_dataset('body' + str(key) + '/body/k',data=data[key].k)
            k.attrs['units'] = ''
            k.attrs['description'] = 'Hydrostatic stiffness matrix'  
            
            
            cg = f.create_dataset('body' + str(key) + '/body/cg',data=data[key].cg)
            cg.attrs['units'] = 'm'
            cg.attrs['description'] = 'Center of gravity'  

            cb = f.create_dataset('body' + str(key) + '/body/cb',data=data[key].cb)
            cb.attrs['units'] = 'm'
            cb.attrs['description'] = 'Center of buoyancy'  
            
            exMag = f.create_dataset('body' + str(key) + '/hydro/ex/mag',data=data[key].ex.mag)
            exMag.attrs['units'] = ''
            exMag.attrs['description'] = 'Magnitude of excitation force'  
            
            exPhase = f.create_dataset('body' + str(key) + '/hydro/ex/phase',data=data[key].ex.phase)
            exPhase.attrs['units'] = 'rad'
            exPhase.attrs['description'] = 'Phase angle of exctiation force'  
            
            exRe = f.create_dataset('body' + str(key) + '/hydro/ex/re',data=data[key].ex.re)
            exRe.attrs['units'] = ''
            exRe.attrs['description'] = 'Real component of excitation force'  

            exIm = f.create_dataset('body' + str(key) + '/hydro/ex/im',data=data[key].ex.im)
            exIm.attrs['units'] = ''
            exIm.attrs['description'] = 'Imaginary component of excitation force'  

            # Write added mass information                
            amInf = f.create_dataset('body' + str(key) + '/hydro/am/inf',data=data[key].am.infFreq)
            amInf.attrs['units for translational degrees of freedom'] = 'kg'
            amInf.attrs['description'] = 'Infinite frequency added mass'
            
            am = f.create_dataset('body' + str(key) + '/hydro/am/all',data=data[key].am.all)
            am.attrs['units for translational degrees of freedom'] = 'kg'                
            am.attrs['units for rotational degrees of freedom'] = 'kg-m^2'
            am.attrs['description'] = 'Added mass. Frequency is the thrid dimension of the data structure.'
            
            for m in xrange(np.shape(data[key].am.all)[0]):
            
                for n in xrange(np.shape(data[key].am.all)[1]):

                    amComp = f.create_dataset('body' + str(key) + '/hydro/am/comps/' + str(m) + ',' + str(n),data=data[key].am.all[m,n,:])
                    amComp.attrs['units'] = ''
                    amComp.attrs['description'] = 'Added mass components as a function of frequency'

                    radComp = f.create_dataset('body' + str(key) + '/hydro/rd/comps/' + str(m) + ',' + str(n),data=data[key].rd.all[m,n,:])
                    radComp.attrs['units'] = ''
                    radComp.attrs['description'] = 'Radiation damping components as a function of frequency'
            
            rad = f.create_dataset('body' + str(key) + '/hydro/rd/all',data=data[key].rd.all)
            rad.attrs['units'] = ''
            rad.attrs['description'] = 'Radiation damping. Frequency is the thrid dimension of the data structure.'
            
            wDepth = f.create_dataset('body' + str(key) + '/sim/wDepth',data=data[key].waterDepth)
            wDepth.attrs['units'] = 'm'
            wDepth.attrs['description'] = 'Water depth'

            waveHead = f.create_dataset('body' + str(key) + '/sim/wDir',data=data[key].waveDir)
            waveHead.attrs['units'] = 'rad'
            waveHead.attrs['description'] = 'Wave direction'
            
            vol = f.create_dataset('body' + str(key) + '/body/dispVol',data=data[key].volDisp)
            vol.attrs['units'] = 'm^3'
            vol.attrs['description'] = 'Displaced volume'

            
            g = f.create_dataset('body' + str(key) + '/sim/g',data=data[key].g)
            g.attrs['units'] = 'm/s^2'
            g.attrs['description'] = 'Gravitational acceleration'
            
            rho = f.create_dataset('body' + str(key) + '/sim/rho',data=data[key].rho)
            rho.attrs['units'] = 'kg/m^3'
            rho.attrs['description'] = 'Water density'
            
            name = f.create_dataset('body' + str(key) + '/sim/name',data=data[key].name)
            name.attrs['description'] = 'Body name'
            
        print 'Wrote HDF5 data to ' + outFile

def writeWecSimHydroData(data,outFile):

    for i in range(np.size(data.keys())):

        import scipy.io as sio
        curData = data[i]
        out = {}
        out['waterDepth'] = curData.waterDepth
        out['waveHeading'] = curData.waveDir
        out['vol'] = curData.volDisp
        out['cg'] = curData.cg
        out['period'] = curData.T[::-1]
        out['linearHyroRestCoef'] = curData.k
        out['fAddedMassZero'] = curData.am.infFreq
        out['fAddedMass'] = curData.am.all[:,:,::-1]
        out['fDamping'] = curData.rd.all[:,:,::-1]
        out['fExtRe'] = curData.ex.re[::-1,:].transpose()
        out['fExtIm'] = curData.ex.im[::-1,:].transpose()
        out['fExtMag'] = curData.ex.mag[::-1,:].transpose()
        out['fExtPhase'] = curData.ex.phase[::-1,:].transpose()
            
        outFileName = outFile[0:-4] + '-body' + str(i) +'.mat'
        sio.savemat(outFileName,out)

        print 'Wrote MATLAB output for WEC-Sim to ' + outFileName

def generateFileNames(outFile):
    '''
    Function to generate filenames needed by hydroData module

    Inputs:
    outFile -- Name of hydrodynamic data file

    Outputs:
    files -- a dictionary of file generateFileNames
    '''

    (path,file) = os.path.split(outFile)
 

    files = {}
    files['out'] = os.path.join(path,file)
    files['hdf5'] = os.path.join(path,file[0:-4] + '.h5')
    files['pickle'] = os.path.join(path,file[0:-4] + '.p')
    files['wecSim'] = os.path.join(path,file[0:-4] + '.mat')

    return files
