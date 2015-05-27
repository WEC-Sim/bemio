# Copyright 2014 the National Renewable Energy Laboratory

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This class contains a structure to store hydrodynamic data from WAMTI,
# AQWA, Nemoh, or another code that calculates hydrodynamic coefficinets
# and excitation forces

# Author: Michael Lawson


import numpy as np

import os

import matplotlib.pyplot as plt

import pickle

import h5py

from scipy import interpolate

from scipy.linalg import hankel, expm

from progressbar import ProgressBar, Bar, Percentage

class Raw(object):
    '''
    Empty class to store raw data from the bemio.io objects
    '''
    def __init__(self):
        pass


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
        self.t = np.array([])
        self.w = np.array([])
        self.L = np.array([])
        self.K = np.array([])

class StateSpaceCoefficients(object):
    '''
    Data class that contains hydrodynamic excitation coefficinets.
    Variables:
    A, B, C, D -- State Space Coefficients
    6 x 6*num_bodies
    '''
    def __init__(self):

        self.A             = np.array([])
        self.B             = np.array([])
        self.C             = np.array([])
        self.D             = np.array([])
        self.irk_bss       = np.array([])
        self.r2t           = np.array([])
        self.rad_conv      = np.array([]) 
        self.it            = np.array([]) 

class HydrodynamicData(object):
    ''''
    Sturcuture for storing data from WAMIT, AQWA and Nemoh.

    Variables:
    rho -- density
    g -- gravity
    files -- Python dictionary of files associated with the
    simulation
    num_bodies -- Total number of bodies in simulation
    body_num -- Body number of the rigid body in the simulation
    cg -- Center of gravity
    cb -- Center of buoyancy
    disp_vol -- Volume displacement
    T -- Wave periods of simulations (e.g. [1, 2, 3, 4, 5,])
    w -- Wave freqencies of simulations 
    am -- Added mass coefficients of HydrodynamicCoefficients type
    rd -- Radiation damping coefficients of HydrodynamicCoefficients
    type
    wp_area -- Water plane area          
    buoy_force -- Buoyancy force at equilibrium
    k -- Hydrostatic stiffness matrix
    ex -- Excitation coeffs of HydrodynamicExcitation type
    water_depth -- Water depth
    wave_dir -- Wave direction 
    name -- Name of the body in the simulation
    '''

    def __init__(self):
        # Default values
        self.rho            = 1000.
        self.g              = 9.81
        self.wave_dir       = 0.
        self.num_bodies     = 0     
                     
        # np.array([])     
        self.cg             = 'not_defined'                          
        self.cb             = 'not_defined'                           
        self.k              = 'not_defined'                           
        self.T              = 'not_defined'                                                       
        self.w              = 'not_defined'
        self.wave_heading   = 'not_defined'                          
        
        # np.floats()
        self.wp_area        = 'not_defined'                             
        self.buoy_force     = 'not_defined'     
        self.disp_vol       = 'not_defined'                         
        self.water_depth    = 'not_defined'                           
        self.body_num       = 'not_defined'                    
        
        # strings
        self.name            = 'not_defined'
        self.bem_code        = 'not_defined'
        self.bem_raw_data    = 'not_defined'

        # objects
        self.am              = HydrodynamicCoefficients()    
        self.rd              = HydrodynamicCoefficients()  
        self.ex              = HydrodynamicExcitation()  
        self.rao             = HydrodynamicExcitation() # same format as hydrodynamic excitation 
        self.ssy             = HydrodynamicExcitation() # same format as hydrodynamic excitation 
        self.irf             = IRF()
        self.ss              = StateSpaceCoefficients()

    def __repr__(self):
        out_string = 'Body name: ' + str(self.name) + \
            '\n    Body number: ' + str(self.body_num) +\
            '\n    Total number of bodies: ' + str(self.num_bodies) + \
            '\n    Displaced volume (m^3): ' + str(self.disp_vol) + \
            '\n    Center of gravity (m): ' + str(self.cg) + \
            '\n    Center of buoyancy (m): ' + str(self.cb)
        return out_string


    def calc_irf(self, t_end=100, n_t = 1001, n_w=1001):
        '''
        Calculate the impulse response functions. See WAMITv7 manual section 13-8

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
        pbar = ProgressBar(widgets=['Calculating the impulse response function for ' + self.name + ':',Percentage(), Bar()], maxval=np.size(self.irf.t)*shape_rd[0]*shape_rd[1]).start()
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

    def calc_ss(self, max_order=10, r2_thresh=0.95 ):
        '''
        Function to calculate state space coefficients
        
        Inputs:
        Kr       - impulse response function
        ss_max    - maximum order of the state space realization
        R2Thresh - R2 threshold that must be met either by the R2 value for K_{r}
        dt       - time step used for the sampling frequency of the impulse response function

        Outputs:
        Ass - time-invariant state matrix
        Bss - time-invariant input matrix
        Css - time-invariant output matrix
        Dss - time-invariant feedthrough matrix
        k_ss_est - Impusle response function as cacluated from state space approximation
        status - status of the realization, 0 - zero hydrodynamic coefficients
        1 - state space realization meets R2 threshold
        2 - state space realization does not
        meet R2 threshold and at ss_max limit
               
        [Ass,Bss,Css,Dss,Krest,status]       
        SS_TD(bodyTemp.hydroForce.irkb(:,ii,jj),simu.ss_max,simu.R2Thresh,simu.dt)
        '''
        dt                  = self.irf.t[2]-self.irf.t[1]
        numFreq             = np.size(self.irf.t) 
        self.ss.irk_bss         = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[0],numFreq] )
        r2bt                = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[0],numFreq] )
        # rd_est    = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[0],numFreq] )
        # am_est  = np.zeros( [ np.shape(self.am.inf)[0],np.shape(self.am.inf)[0],numFreq] )
        k_ss_est               = np.zeros( numFreq )
        self.ss.A = np.zeros([6,np.shape(self.am.inf)[1],max_order,max_order])
        self.ss.B = np.zeros([6,np.shape(self.am.inf)[1],max_order,1])
        self.ss.C = np.zeros([6,np.shape(self.am.inf)[1],1,max_order])
        self.ss.D = np.zeros([6,np.shape(self.am.inf)[1],1])
        self.ss.irk_bss   = np.zeros([6,np.shape(self.am.inf)[1],numFreq])
        self.ss.rad_conv= np.zeros([6,np.shape(self.am.inf)[1]])
        self.ss.it     = np.zeros([6,np.shape(self.am.inf)[1]])
        self.ss.r2t      = np.zeros([6,np.shape(self.am.inf)[1]])
        
        pbar = ProgressBar(widgets=['Calculating state space coefficients for ' + self.name + ':',Percentage(), Bar()], maxval=np.shape(self.am.inf)[0]*np.shape(self.am.inf)[1]).start()
        count = 0
        for i in xrange(np.shape(self.am.inf)[0]):

            for j in xrange(np.shape(self.am.inf)[1]):

                r2bt = np.linalg.norm(self.irf.K[i,j,:]-self.irf.K.mean(axis=2)[i,j])
                
                ss = 2 #Initial state space order
                while True:
                    
                    #Perform Hankel Singular Value Decomposition
                    y=dt*self.irf.K[i,j,:]                    
                    h=hankel(y[1::])
                    u,svh,v=np.linalg.svd(h)
                    
                    u1 = u[0:numFreq-2,0:ss]
                    v1 = v.T[0:numFreq-2,0:ss]
                    u2 = u[1:numFreq-1,0:ss]
                    sqs = np.sqrt(svh[0:ss].reshape(ss,1))
                    invss = 1/sqs
                    ubar = np.dot(u1.T,u2)

                    a = ubar*np.dot(invss,sqs.T)
                    b = v1[0,:].reshape(ss,1)*sqs
                    c = u1[0,:].reshape(1,ss)*sqs.T
                    d = y[0]        

                    CoeA = dt/2
                    CoeB = 1
                    CoeC = -CoeA
                    CoeD = 1

                    iidd = np.linalg.inv(CoeA*np.eye(ss)-CoeC*a)               #(T/2*I + T/2*A)^{-1}         = 2/T(I + A)^{-1}
                    
                    ac = np.dot(CoeB*a-CoeD*np.eye(ss),iidd)                            #(A-I)2/T(I + A)^{-1}         = 2/T(A-I)(I + A)^{-1}
                    bc = (CoeA*CoeB-CoeC*CoeD)*np.dot(iidd,b)                         #(T/2+T/2)*2/T(I + A)^{-1}B   = 2(I + A)^{-1}B
                    cc = np.dot(c,iidd)                                                #C * 2/T(I + A)^{-1}          = 2/T(I + A)^{-1}
                    dc = d + CoeC*np.dot(np.dot(c,iidd),b)                                     #D - T/2C (2/T(I + A)^{-1})B  = D - C(I + A)^{-1})B

                    for jj in xrange(numFreq):

                        k_ss_est[jj] = np.dot(np.dot(cc,expm(ac*dt*jj)),bc)                         #Calculate impulse response function from state space approximation
  
                    R2TT = np.linalg.norm(self.irf.K[i,j,:]-k_ss_est)                #Calculate 2 norm of the difference between know and estimated values impulse response function
                    R2T = 1 - np.square(R2TT/r2bt)                                    #Calculate the R2 value for impulse response function

                    if R2T >= r2_thresh:                                   #Check to see if threshold for the impulse response is meet
                    
                        status = 1                                             #%Set status
                        break
                    
                    if ss == max_order:                                            #Check to see if limit on the state space order has been reached
                    
                        status = 2                                             #%Set status
                        break
                    
                    ss=ss+1                                                    #Increase state space order
                                        
                self.ss.A[i,j,0:np.shape(ac)[0],0:np.shape(ac)[0]]  = ac
                self.ss.B[i,j,0:np.shape(bc)[0],0                ]  = bc[:,0]
                self.ss.C[i,j,0                ,0:np.shape(cc)[1]]  = cc[0,:]
                self.ss.D[i,j]                                      = dc
                self.ss.irk_bss[i,j,:]  = k_ss_est
                self.ss.rad_conv[i,j] = status
                self.ss.r2t[i,j] = R2T
                self.ss.it[i,j] = ss
                count += 1
                pbar.update(count)

        pbar.finish()
        
    def plot_irf(self,components):
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
        

    def plot_am_rd(self,components):
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

    def plot_excitation(self,components):
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

            ax[0].plot(w,mag,'x-',label='Component (' + str(m+1) + ')')
            ax[1].plot(w,phase,'x-',label='Component (' + str(m+1) + ')')
            ax[2].plot(w,re,'x-',label='Component (' + str(m+1) + ')')
            ax[3].plot(w,im,'x-',label='Component (' + str(m+1) + ')')

            ax[0].legend(loc=0)


def write_pickle(data_obj,out_file=None):
    '''
    Writes hydrodynamic data to a pickle file.
    
    Inputs:
    data -- dictionary that contains HydrodynamicData objects for each body in the simulation
    out_file -- name of the pickle file
    '''
    if out_file is None:
        out_file = data_obj.files['pickle']

    pickle.dump(data_obj,open(out_file,'wb'))
    
    print 'Wrote pickle data to ' + out_file


def write_hdf5(data_obj,out_file=None):
    '''
    Writes hydrodynamic data to a HDF5 file structure.
    
    Inputs:
    data -- dictionary that contains HydrodynamicData objects for each body in the simulation
    outFile -- name of the hdf5 file

    Outputs: None
    '''

    if out_file is None:
        out_file = data_obj.files['hdf5']


    with h5py.File(out_file, "w") as f:       

        for key, key in enumerate(data_obj.data.keys()):

            # Body properities
            cg = f.create_dataset('body' + str(key+1) + '/properties/cg',data=data_obj.data[key].cg)
            cg.attrs['units'] = 'm'
            cg.attrs['description'] = 'Center of gravity'  

            cb = f.create_dataset('body' + str(key+1) + '/properties/cb',data=data_obj.data[key].cb)
            cb.attrs['units'] = 'm'
            cb.attrs['description'] = 'Center of buoyancy' 

            vol = f.create_dataset('body' + str(key+1) + '/properties/disp_vol',data=data_obj.data[key].disp_vol)
            vol.attrs['units'] = 'm^3'
            vol.attrs['description'] = 'Displaced volume'

            name = f.create_dataset('body' + str(key+1) + '/properties/name',data=data_obj.data[key].name)
            name.attrs['description'] = 'Name of rigid body'

            num = f.create_dataset('body' + str(key+1) + '/properties/body_number',data=data_obj.data[key].body_num)
            num.attrs['description'] = 'Number of rigid body from the BEM simulation'
            
            # Hydro coeffs
            try:

                irfK = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/K',data=data_obj.data[key].irf.K)
                irfK.attrs['units'] = ''
                irfK.attrs['description'] = 'Impulse response function' 

                irfT = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/t',data=data_obj.data[key].irf.t)
                irfT.attrs['units'] = 'seconds'
                irfT.attrs['description'] = 'Time vector for the impulse response function' 

                irfW = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/w',data=data_obj.data[key].irf.w)
                irfW.attrs['units'] = 'seconds'
                irfW.attrs['description'] = 'Interpolated frequencies used to compute the impulse response function' 

                irfL = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/L',data=data_obj.data[key].irf.L)
                irfL.attrs['units'] = ''
                irfL.attrs['description'] = 'Time derivatitive of the impulse response functiuon' 

                for m in xrange(np.shape(data_obj.data[key].am.all)[0]):
            
                    for n in xrange(np.shape(data_obj.data[key].am.all)[1]):

                        irfLComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/comps/L/comp_' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].irf.L[m,n,:])
                        irfLComp.attrs['units'] = ''
                        irfLComp.attrs['description'] = 'Components of the IRF'

                        irfKComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/comps/K/comp_' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].irf.K[m,n,:])
                        irfKComp.attrs['units'] = ''
                        irfKComp.attrs['description'] = 'Components of the ddt(IRF): K'
            except:

                print 'IRF functions for ' + data_obj.data[key].name + ' were not written because they were not calculated. Use the calc_irf function to calculate the IRF.'

            try:
    
                ssRadfA = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/A/all',data=data_obj.data[key].ss.A)
                ssRadfA.attrs['units'] = ''
                ssRadfA.attrs['description'] = 'State Space A Coefficient'
                
                ssRadfB = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/B/all',data=data_obj.data[key].ss.B)
                ssRadfB.attrs['units'] = ''
                ssRadfB.attrs['description'] = 'State Space B Coefficient'

                ssRadfC = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/C/all',data=data_obj.data[key].ss.C)
                ssRadfC.attrs['units'] = ''
                ssRadfC.attrs['description'] = 'State Space C Coefficient'

                ssRadfD = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/D/all',data=data_obj.data[key].ss.D)
                ssRadfD.attrs['units'] = ''
                ssRadfD.attrs['description'] = 'State Space D Coefficient'

                r2t = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/r2t',data=data_obj.data[key].ss.r2t)
                r2t.attrs['units'] = ''
                r2t.attrs['description'] = 'State space curve fitting R**2 value'

                it = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/it',data=data_obj.data[key].ss.it)
                it.attrs['units'] = ''
                it.attrs['description'] = 'Order of state space reailization'

                for m in xrange(np.shape(data_obj.data[key].am.all)[0]):
            
                    for n in xrange(np.shape(data_obj.data[key].am.all)[1]):

                        ss_A = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/A/comps/comp_' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].ss.A[m,n,:,:])
                        ss_A.attrs['units'] = ''
                        ss_A.attrs['description'] = 'Components of the State Space A Coefficient'

                        ss_B = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/B/comps/comp_' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].ss.B[m,n,:,:])
                        ss_B.attrs['units'] = ''
                        ss_B.attrs['description'] = 'Components of the State Space B Coefficient'

                        ss_C = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/C/comps/comp_' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].ss.C[m,n,:,:])
                        ss_C.attrs['units'] = ''
                        ss_C.attrs['description'] = 'Components of the State Space C Coefficient'

                        ss_D = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/D/comps/comp_' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].ss.D[m,n])
                        ss_D.attrs['units'] = ''
                        ss_D.attrs['description'] = 'Components of the State Space C Coefficient'
    
            except:

                print 'State Space Coefficients for ' + data_obj.data[key].name + ' were not written because they were not calculated. Use the calcSS function to calculate the State Space Coefficients.'

            k = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/linear_restoring_stiffness',data=data_obj.data[key].k)
            k.attrs['units'] = ''
            k.attrs['description'] = 'Hydrostatic stiffness matrix'  

            exMag = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/mag',data=data_obj.data[key].ex.mag)
            exMag.attrs['units'] = ''
            exMag.attrs['description'] = 'Magnitude of excitation force'  
            
            exPhase = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/phase',data=data_obj.data[key].ex.phase)
            exPhase.attrs['units'] = 'rad'
            exPhase.attrs['description'] = 'Phase angle of exctiation force'  
            
            exRe = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/re',data=data_obj.data[key].ex.re)
            exRe.attrs['units'] = ''
            exRe.attrs['description'] = 'Real component of excitation force'  

            exIm = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/im',data=data_obj.data[key].ex.im)
            exIm.attrs['units'] = ''
            exIm.attrs['description'] = 'Imaginary component of excitation force'  

            # Write added mass information                
            amInf = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/inf_freq',data=data_obj.data[key].am.inf)
            amInf.attrs['units for translational degrees of freedom'] = 'kg'
            amInf.attrs['description'] = 'Infinite frequency added mass'
            
            am = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/all',data=data_obj.data[key].am.all)
            am.attrs['units for translational degrees of freedom'] = 'kg'                
            am.attrs['units for rotational degrees of freedom'] = 'kg-m^2'
            am.attrs['description'] = 'Added mass. Frequency is the thrid dimension of the data structure.'
            
            for m in xrange(np.shape(data_obj.data[key].am.all)[0]):
            
                for n in xrange(np.shape(data_obj.data[key].am.all)[1]):

                    amComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/comps/comp_' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].am.all[m,n,:])
                    amComp.attrs['units'] = ''
                    amComp.attrs['description'] = 'Added mass components as a function of frequency'

                    radComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/comps/' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].rd.all[m,n,:])
                    radComp.attrs['units'] = ''
                    radComp.attrs['description'] = 'Radiation damping components as a function of frequency'
            
            rad = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/all',data=data_obj.data[key].rd.all)
            rad.attrs['units'] = ''
            rad.attrs['description'] = 'Radiation damping. Frequency is the thrid dimension of the data structure.'

        # Simulation parameters
        g = f.create_dataset('simulation_parameters/g',data=data_obj.data[key].g)
        g.attrs['units'] = 'm/s^2'
        g.attrs['description'] = 'Gravitational acceleration'
        
        rho = f.create_dataset('simulation_parameters/rho',data=data_obj.data[key].rho)
        rho.attrs['units'] = 'kg/m^3'
        rho.attrs['description'] = 'Water density'

        T = f.create_dataset('simulation_parameters/T',data=data_obj.data[key].T)
        T.attrs['units'] = 's'
        T.attrs['description'] = 'Wave periods'
        
        w = f.create_dataset('simulation_parameters/w',data=data_obj.data[key].w)
        w.attrs['units'] = 'rad/s'                
        w.attrs['description'] = 'Wave frequencies'

        water_depth = f.create_dataset('simulation_parameters/water_depth',data=data_obj.data[key].water_depth)
        water_depth.attrs['units'] = 'm'
        water_depth.attrs['description'] = 'Water depth'

        wave_dir = f.create_dataset('simulation_parameters/wave_dir',data=data_obj.data[key].wave_dir)
        wave_dir.attrs['units'] = 'rad'
        wave_dir.attrs['description'] = 'Wave direction'

        rawOut = f.create_dataset('bem_data/output_file',data=data_obj.data[key].bem_raw_data)
        rawOut.attrs['description'] = 'Raw output from BEM code'

        code = f.create_dataset('simulation_parameters/bem_code',data=data_obj.data[key].bem_code)
        code.attrs['description'] = 'BEM code'

        print 'Wrote HDF5 data to ' + out_file


def generate_file_names(out_file):
    '''
    Function to generate filenames needed by hydroData module

    Inputs:
    outFile -- Name of hydrodynamic data file

    Outputs:
    files -- a dictionary of file generate_file_names
    '''
    out_file = os.path.abspath(out_file)
    (path,file) = os.path.split(out_file)
 
    files = {}
    files['out'] = os.path.join(path,file)
    files['hdf5'] = os.path.join(path,file[0:-4] + '.h5')
    files['pickle'] = os.path.join(path,file[0:-4] + '.p')

    return files
