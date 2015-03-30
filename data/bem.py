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

import pickle

from scipy import interpolate

from progressbar import ProgressBar, Bar, Percentage

from scipy.linalg import hankel, expm

#from scipy.signal import ss2tf


class ViscousDamping(object):
    '''
    This data contains data that defines the viscous damping 
    properities of rigid bodies

    Variables:
    cd -- drag coefficients for the 6 degrees of freedom
    a -- characteristic area for the 6 degrees of freedom
    '''
    def __init__(self):

        self.cd = np.zeros(6)
        self.a  = np.zeros(6)


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

class StateSpaceCoefficients(object):
    '''
    Data class that contains hydrodynamic excitation coefficinets.
    Variables:
    A, B, C, D -- State Space Coefficients
    6 x 6
 
   '''

    def __init__(self):

        self.A             = np.array([])
        self.B             = np.array([])
        self.C             = np.array([])
        self.D             = np.array([])


class IRF(object):
    '''
    Object that contains the IRF data

    Variables:
    dt -- timestep for the IRF calculation
    t_end -- end time for the IRF calculation
    t_series -- time series for the IRF calculation
    L -- impulse response function
    K -- time derivavitative of the impulse response function
    '''

    def __init__(self):
        self.t = None
        self.w = None
        self.L = None
        self.K = None


class HydrodynamicData(object):
    ''''
    Sturcuture for storing data from WAMIT, AQWA and Nemoh.

    Variables:
    rho -- density
    g -- gravity
    files -- Python dictionary of files associated with the
    simulation
    nBodies -- Total number of bodies in simulation
    bodyN -- Body number of the rigid body in the simulation
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
        self.waveDir        = 0.
        self.nBodies        = 0     
                          
        self.files          = {}
        self.cg             = {}                            
        self.cb             = {}                            
        self.volDisp        = {}                            
        self.T              = {}                            
        self.w              = {}                            
        self.wpArea         = {}                            
        self.buoyForce      = {}                            
        self.k              = {}                            
           
        self.waterDepth     = None                          
        self.bodyN          = None                      
        self.name           = None

        self.am             = HydrodynamicCoefficients()    
        self.rd             = HydrodynamicCoefficients()  
        self.ex             = HydrodynamicExcitation()   
        self.vDamping       = ViscousDamping()
        self.irf            = IRF()
        self.ssRadf         = StateSpaceCoefficients()    
        self.ssMax          = 10
        self.R2Thresh       = 0.95
        self.irkbss         = np.array([])
        self.R2T            = np.array([])
        self.ssRadconv      = np.array([])
        self.ssIt           = np.array([])
#        self.fDampingest    = np.array([])
#        self.fAddedMassest  = np.array([])
    
    def calcIRF(self, t_end=100, n_t = 10001, n_w=1001):
        '''
        Calculate the IRF. See WAMITv7 manual section 13-8

        Inputs:

        Outputs:
        This function populates the irf variable
        '''

        self.irf.t = np.linspace(0,t_end,n_t)
        self.irf.w = np.linspace(np.min(self.w),np.max(self.w),n_w)

        self.irf.L = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[1],np.size(self.irf.t) ] )
        self.irf.K = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[1],np.size(self.irf.t) ] )

        rd_interp = np.zeros( [ np.shape(self.rd.all)[0], np.shape(self.rd.all)[1], np.size(self.irf.w) ])

        shape_rd = np.shape(self.rd.all)

        # Interpolate the radiation damping matrix for the IRF calculation
        flip = False

        if self.w[0] > self.w[1]:

            wTmp = np.flipud(self.w)
            flip = True

        else:

            wTmp = self.w

        for i in xrange(shape_rd[0]):

            for j in xrange(shape_rd[1]):

                if flip is True:

                    rdTmp = np.flipud(self.rd.all[i,j,:])

                else:
                    rdTmp = self.rd.all[i,j,:]

                f = interpolate.interp1d(x=wTmp, y=rdTmp)
                rd_interp[i,j,:] = f(self.irf.w) 

        # Calculate the IRF
        pbar = ProgressBar(widgets=['Calculating IRF for ' + self.name + ':',Percentage(), Bar()], maxval=np.size(self.irf.t)*shape_rd[0]*shape_rd[1]).start()
        count = 1
        for t_ind, t in enumerate(self.irf.t):

            for i in xrange(shape_rd[0]):

                for j in xrange(shape_rd[1]):
                    # Radiation damping calculation method
                    tmpL = 2./np.pi*rd_interp[i,j,:]*np.sin(self.irf.w*t)
                    tmpK = 2./np.pi*rd_interp[i,j,:]*np.cos(self.irf.w*t)
                    self.irf.K[i,j,t_ind] = np.trapz(y=tmpK,x=self.irf.w)
                    self.irf.L[i,j,t_ind] = np.trapz(y=tmpL,x=self.irf.w)
                    pbar.update(count)
                    count += 1
 
        pbar.finish()

    def calcSS(self):
        '''
        Function to calculate state space coefficients
        
        Inputs:
        Kr       - impulse response function
        ssMax    - maximum order of the state space realization
        R2Thresh - R2 threshold that must be met either by the R2 value for K_{r}
        dt       - time step used for the sampling frequency of the impulse response function

        Outputs:
        Ass - time-invariant state matrix
        Bss - time-invariant input matrix
        Css - time-invariant output matrix
        Dss - time-invariant feedthrough matrix
        Krlti - Impusle response function as cacluated from state space approximation
        status - status of the realization, 0 - zero hydrodynamic coefficients
                                            1 - state space realization meets R2 threshold
                                            2 - state space realization does not
                                                meet R2 threshold and at ssMax limit
        
        [Ass,Bss,Css,Dss,Krest,status]
        
        SS_TD(bodyTemp.hydroForce.irkb(:,ii,jj),simu.ssMax,simu.R2Thresh,simu.dt)

        '''
        dt                  = self.irf.t[2]-self.irf.t[1]
        numFreq             = np.size(self.irf.t) 
        self.irkbss         = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[0],numFreq] )
        R2BT                = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[0],numFreq] )
        self.fDampingest    = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[0],numFreq] )
        self.fAddedMassest  = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[0],numFreq] )
        Krlti               = np.zeros( numFreq )
        self.ssRadf.A = np.zeros([6,np.shape(self.am.inf)[1],self.ssMax,self.ssMax])
        self.ssRadf.B = np.zeros([6,np.shape(self.am.inf)[1],self.ssMax,1])
        self.ssRadf.C = np.zeros([6,np.shape(self.am.inf)[1],1,self.ssMax])
        self.ssRadf.D = np.zeros([6,np.shape(self.am.inf)[1],1])
        self.irkbss   = np.zeros([6,np.shape(self.am.inf)[1],numFreq])
        self.ssRadconv= np.zeros([6,np.shape(self.am.inf)[1]])
        self.ssIt     = np.zeros([6,np.shape(self.am.inf)[1]])
        self.R2T      = np.zeros([6,np.shape(self.am.inf)[1]])
        
        for i in xrange(np.shape(self.am.inf)[0]):

            for j in xrange(np.shape(self.am.inf)[1]):

                R2BT = np.linalg.norm(self.irf.K[i,j,:]-self.irf.K.mean(axis=2)[i,j])
                
                ss = 2 #Initial state space order
                while True:
                    
                    #Perform Hankel Singular Value Decomposition
                    y=dt*self.irf.K[i,j,:]                    
                    h=hankel(y[1::])
                    u,svh,v=np.linalg.svd(h)
                    
#                    print np.shape(h),np.shape(u),np.shape(svh),np.shape(v)
                    u1 = u[0:numFreq-2,0:ss]
                    v1 = v.T[0:numFreq-2,0:ss]
                    u2 = u[1:numFreq-1,0:ss]
                    sqs = np.sqrt(svh[0:ss].reshape(ss,1))
                    invss = 1/sqs
                    ubar = np.dot(u1.T,u2)

#                    print np.shape(ubar),np.shape(invss),np.shape(sqs)

                    a = ubar*np.dot(invss,sqs.T)
                    b = v1[0,:].reshape(ss,1)*sqs
                    c = u1[0,:].reshape(1,ss)*sqs.T
                    d = y[0]        
#                    print np.shape(a),np.shape(b),np.shape(c),np.shape(d)

                    CoeA = dt/2
                    CoeB = 1
                    CoeC = -CoeA
                    CoeD = 1

                    iidd = np.linalg.inv(CoeA*np.eye(ss)-CoeC*a)               #(T/2*I + T/2*A)^{-1}         = 2/T(I + A)^{-1}
                    
                    ac = np.dot(CoeB*a-CoeD*np.eye(ss),iidd)                            #(A-I)2/T(I + A)^{-1}         = 2/T(A-I)(I + A)^{-1}
                    bc = (CoeA*CoeB-CoeC*CoeD)*np.dot(iidd,b)                         #(T/2+T/2)*2/T(I + A)^{-1}B   = 2(I + A)^{-1}B
                    cc = np.dot(c,iidd)                                                #C * 2/T(I + A)^{-1}          = 2/T(I + A)^{-1}
                    dc = d + CoeC*np.dot(np.dot(c,iidd),b)                                     #D - T/2C (2/T(I + A)^{-1})B  = D - C(I + A)^{-1})B
#                    print np.shape(iidd),np.shape(ac),np.shape(bc),np.shape(cc),np.shape(dc)

                    for jj in xrange(numFreq):                 
                        Krlti[jj] = np.dot(np.dot(cc,expm(ac*dt*jj)),bc)                         #Calculate impulse response function from state space approximation
  
                    R2TT = np.linalg.norm(self.irf.K[i,j,:]-Krlti)                #Calculate 2 norm of the difference between know and estimated values impulse response function

                    R2T = 1 - np.square(R2TT/R2BT)                                    #Calculate the R2 value for impulse response function

                    if R2T >= self.R2Thresh:                                   #Check to see if threshold for the impulse response is meet
                        status = 1                                             #%Set status
                        break
                    if ss == self.ssMax:                                            #Check to see if limit on the state space order has been reached
                        status = 2                                             #%Set status
                        break
                    
                    ss=ss+1                                                    #Increase state space order
                                        
                print i,j,status,np.shape(self.ssRadf.A[i,j,:,:]),np.shape(ac)
                
                self.ssRadf.A[i,j,0:np.shape(ac)[0],0:np.shape(ac)[0]]  = ac
                self.ssRadf.B[i,j,0:np.shape(bc)[0],0                ]  = bc[:,0]
                self.ssRadf.C[i,j,0                ,0:np.shape(cc)[1]]  = cc[0,:]
                self.ssRadf.D[i,j]                                      = dc
                self.irkbss[i,j,:]  = Krlti
                self.ssRadconv[i,j] = status
                self.R2T[i,j] = R2T
                self.ssIt[i,j] = ss
                
        print self.R2T
        print self.ssRadconv

    def plotIRF(self,components):
        '''
        Function to plot the IRF

        Inputs:
        components -- A list of components to plot. E.g [[0,0],[1,1],[2,2]]
        
        Outputs:
        None -- A plot is displayed. The plt.show() command may need to be used
        depending on your python env settings
        '''  
        
        f, ax = plt.subplots(np.shape(components)[0], sharex=True, figsize=(8,10))
                
        # Plot added mass and damping
        for i,comp in enumerate(components):
            
            x = comp[0]
            y = comp[1]
            t = self.irf.t
            L = self.irf.L[x,y,:]
            K = self.irf.K[x,y,:]

            ax[i].set_ylabel('comp ' + str(x) + ',' + str(y))

            ax[i].plot(t,L,label='L')
            ax[i].plot(t,K,label='K ddt(L)')
                  
        ax[0].set_title('IRF for ' + str(self.name))
        ax[0].legend()
        ax[i].set_xlabel('Time (s)')
        

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
        ax[0].set_title('Hydrodynamic coefficients for ' + str(self.name))    
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
        ax[0].set_title('Excitation force for ' + str(self.name))    

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

            # Body properities
            cg = f.create_dataset('body' + str(key) + '/properties/cg',data=data[key].cg)
            cg.attrs['units'] = 'm'
            cg.attrs['description'] = 'Center of gravity'  

            cb = f.create_dataset('body' + str(key) + '/properties/cb',data=data[key].cb)
            cb.attrs['units'] = 'm'
            cb.attrs['description'] = 'Center of buoyancy' 

            vol = f.create_dataset('body' + str(key) + '/properties/dispVol',data=data[key].volDisp)
            vol.attrs['units'] = 'm^3'
            vol.attrs['description'] = 'Displaced volume'

            name = f.create_dataset('body' + str(key) + '/properties/name',data=data[key].name)
            name.attrs['description'] = 'Name of rigid body'

            num = f.create_dataset('body' + str(key) + '/properties/bodyNumber',data=data[key].bodyN)
            num.attrs['description'] = 'Number of rigid body from the BEM simulation'
            
            # Hydro coeffs
            try:
                irfK = f.create_dataset('body' + str(key) + '/hydro_coeffs/irf/K',data=data[key].irf.K)
                irfK.attrs['units'] = ''
                irfK.attrs['description'] = 'Impulse response function' 
    
                irfT = f.create_dataset('body' + str(key) + '/hydro_coeffs/irf/t',data=data[key].irf.t)
                irfT.attrs['units'] = 'seconds'
                irfT.attrs['description'] = 'Time vector for the impulse response function' 
    
                irfW = f.create_dataset('body' + str(key) + '/hydro_coeffs/irf/w',data=data[key].irf.w)
                irfW.attrs['units'] = 'seconds'
                irfW.attrs['description'] = 'Interpolated frequencies used to compute the impulse response function' 
    
                irfL = f.create_dataset('body' + str(key) + '/hydro_coeffs/irf/L',data=data[key].irf.L)
                irfL.attrs['units'] = ''
                irfL.attrs['description'] = 'Time derivatitive of the impulse response functiuon' 
    
                for m in xrange(np.shape(data[key].am.all)[0]):
            
                    for n in xrange(np.shape(data[key].am.all)[1]):
    
                        irfLComp = f.create_dataset('body' + str(key) + '/hydro_coeffs/irf/comps/L/comp_' + str(m) + '_' + str(n),data=data[key].irf.L[m,n,:])
                        irfLComp.attrs['units'] = ''
                        irfLComp.attrs['description'] = 'Components of the IRF'
    
                        irfKComp = f.create_dataset('body' + str(key) + '/hydro_coeffs/irf/comps/K/comp_' + str(m) + '_' + str(n),data=data[key].irf.K[m,n,:])
                        irfKComp.attrs['units'] = ''
                        irfKComp.attrs['description'] = 'Components of the ddt(IRF): K'
        
            except:

                print 'IRF functions for ' + data[key].name + ' were not written because they were not calculated. Use the calcIRF function to calculate the IRF.'

            try:
    
                ssRadfA = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/A',data=data[key].ssRadf.A)
                ssRadfA.attrs['units'] = ''
                ssRadfA.attrs['description'] = 'State Space A Coefficient'
                
                ssRadfB = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/B',data=data[key].ssRadf.B)
                ssRadfB.attrs['units'] = ''
                ssRadfB.attrs['description'] = 'State Space B Coefficient'

                ssRadfC = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/C',data=data[key].ssRadf.C)
                ssRadfC.attrs['units'] = ''
                ssRadfC.attrs['description'] = 'State Space C Coefficient'

                ssRadfD = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/D',data=data[key].ssRadf.D)
                ssRadfD.attrs['units'] = ''
                ssRadfD.attrs['description'] = 'State Space D Coefficient'

                for m in xrange(np.shape(data[key].am.all)[0]):
            
                    for n in xrange(np.shape(data[key].am.all)[1]):
    
                        ssRadfAComp = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/A/comp_' + str(m) + '_' + str(n),data=data[key].ssRadf.A[m,n,:,:])
                        ssRadfAComp.attrs['units'] = ''
                        ssRadfAComp.attrs['description'] = 'Components of the State Space A Coefficient'
    
                        ssRadfBComp = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/B/comp_' + str(m) + '_' + str(n),data=data[key].ssRadf.B[m,n,:,:])
                        ssRadfBComp.attrs['units'] = ''
                        ssRadfBComp.attrs['description'] = 'Components of the State Space B Coefficient'
    
                        ssRadfCComp = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/C/comp_' + str(m) + '_' + str(n),data=data[key].ssRadf.C[m,n,:,:])
                        ssRadfCComp.attrs['units'] = ''
                        ssRadfCComp.attrs['description'] = 'Components of the State Space C Coefficient'
    
                        ssRadfDComp = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/D/comp_' + str(m) + '_' + str(n),data=data[key].ssRadf.D[m,n])
                        ssRadfDComp.attrs['units'] = ''
                        ssRadfDComp.attrs['description'] = 'Components of the State Space C Coefficient'
    
            except:

                print 'State Space Coefficients for ' + data[key].name + ' were not written because they were not calculated. Use the calcSS function to calculate the State Space Coefficients.'

            k = f.create_dataset('body' + str(key) + '/hydro_coeffs/k',data=data[key].k)
            k.attrs['units'] = ''
            k.attrs['description'] = 'Hydrostatic stiffness matrix'  

            exMag = f.create_dataset('body' + str(key) + '/hydro_coeffs/ex/mag',data=data[key].ex.mag)
            exMag.attrs['units'] = ''
            exMag.attrs['description'] = 'Magnitude of excitation force'  
            
            exPhase = f.create_dataset('body' + str(key) + '/hydro_coeffs/ex/phase',data=data[key].ex.phase)
            exPhase.attrs['units'] = 'rad'
            exPhase.attrs['description'] = 'Phase angle of exctiation force'  
            
            exRe = f.create_dataset('body' + str(key) + '/hydro_coeffs/ex/re',data=data[key].ex.re)
            exRe.attrs['units'] = ''
            exRe.attrs['description'] = 'Real component of excitation force'  

            exIm = f.create_dataset('body' + str(key) + '/hydro_coeffs/ex/im',data=data[key].ex.im)
            exIm.attrs['units'] = ''
            exIm.attrs['description'] = 'Imaginary component of excitation force'  

            # Write added mass information                
            amInf = f.create_dataset('body' + str(key) + '/hydro_coeffs/am/inf',data=data[key].am.inf)
            amInf.attrs['units for translational degrees of freedom'] = 'kg'
            amInf.attrs['description'] = 'Infinite frequency added mass'
            
            am = f.create_dataset('body' + str(key) + '/hydro_coeffs/am/all',data=data[key].am.all)
            am.attrs['units for translational degrees of freedom'] = 'kg'                
            am.attrs['units for rotational degrees of freedom'] = 'kg-m^2'
            am.attrs['description'] = 'Added mass. Frequency is the thrid dimension of the data structure.'
            
            for m in xrange(np.shape(data[key].am.all)[0]):
            
                for n in xrange(np.shape(data[key].am.all)[1]):

                    amComp = f.create_dataset('body' + str(key) + '/hydro_coeffs/am/comps/comp_' + str(m) + '_' + str(n),data=data[key].am.all[m,n,:])
                    amComp.attrs['units'] = ''
                    amComp.attrs['description'] = 'Added mass components as a function of frequency'

                    radComp = f.create_dataset('body' + str(key) + '/hydro_coeffs/rd/comps/' + str(m) + '_' + str(n),data=data[key].rd.all[m,n,:])
                    radComp.attrs['units'] = ''
                    radComp.attrs['description'] = 'Radiation damping components as a function of frequency'
            
            rad = f.create_dataset('body' + str(key) + '/hydro_coeffs/rd/all',data=data[key].rd.all)
            rad.attrs['units'] = ''
            rad.attrs['description'] = 'Radiation damping. Frequency is the thrid dimension of the data structure.'

            ssIt = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/ss',data=data[key].ssIt)
            ssIt.attrs['units'] = ''
            ssIt.attrs['description'] = 'State Space Matrix Size'

            ssStatus = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/status',data=data[key].ssRadconv)
            ssStatus.attrs['units'] = ''
            ssStatus.attrs['description'] = 'State Space Convergence Status'

            R2T = f.create_dataset('body' + str(key) + '/hydro_coeffs/stateSpace/R2T',data=data[key].R2T)
            R2T.attrs['units'] = ''
            R2T.attrs['description'] = 'State Space R2T'

#            ssRadfB = f.create_dataset('body' + str(key) + '/hydro_coeffs/ssRadf/B',data=data[key].ssRadf.B)
#            ssRadfB.attrs['units'] = ''
#            ssRadfB.attrs['description'] = 'State Space Coefficients'
#
#            ssRadfC = f.create_dataset('body' + str(key) + '/hydro_coeffs/ssRadf/C',data=data[key].ssRadf.C)
#            ssRadfC.attrs['units'] = ''
#            ssRadfC.attrs['description'] = 'State Space Coefficients'
#
#            ssRadfD = f.create_dataset('body' + str(key) + '/hydro_coeffs/ssRadf/D',data=data[key].ssRadf.D)
#            ssRadfD.attrs['units'] = ''
#            ssRadfD.attrs['description'] = 'State Space Coefficients'

        # Simulation parameters
        g = f.create_dataset('simulation_parameters/g',data=data[key].g)
        g.attrs['units'] = 'm/s^2'
        g.attrs['description'] = 'Gravitational acceleration'
        
        rho = f.create_dataset('simulation_parameters/rho',data=data[key].rho)
        rho.attrs['units'] = 'kg/m^3'
        rho.attrs['description'] = 'Water density'

        T = f.create_dataset('simulation_parameters/T',data=data[key].T)
        T.attrs['units'] = 's'
        T.attrs['description'] = 'Wave periods'
        
        w = f.create_dataset('simulation_parameters/w',data=data[key].w)
        w.attrs['units'] = 'rad/s'                
        w.attrs['description'] = 'Wave frequencies'

        wDepth = f.create_dataset('simulation_parameters/wDepth',data=data[key].waterDepth)
        wDepth.attrs['units'] = 'm'
        wDepth.attrs['description'] = 'Water depth'

        waveHead = f.create_dataset('simulation_parameters/wDir',data=data[key].waveDir)
        waveHead.attrs['units'] = 'rad'
        waveHead.attrs['description'] = 'Wave direction'        

        ssMax = f.create_dataset('simulation_parameters/ssMax',data=data[key].ssMax)
        ssMax.attrs['units'] = ''
        ssMax.attrs['description'] = 'The upper limit on the state space order constructed from realization program'      

        R2Thresh = f.create_dataset('simulation_parameters/R2Thresh',data=data[key].R2Thresh)
        R2Thresh.attrs['units'] = ''
        R2Thresh.attrs['description'] = 'The threshold set on R^2 to stop the realization program'      


            
        print 'Wrote HDF5 data to ' + outFile


def generateFileNames(out_file):
    '''
    Function to generate filenames needed by hydroData module

    Inputs:
    outFile -- Name of hydrodynamic data file

    Outputs:
    files -- a dictionary of file generateFileNames
    '''
    out_file = os.path.abspath(out_file)
    (path,file) = os.path.split(out_file)
 
    files = {}
    files['out'] = os.path.join(path,file)
    files['hdf5'] = os.path.join(path,file[0:-4] + '.h5')
    files['pickle'] = os.path.join(path,file[0:-4] + '.p')

    return files
