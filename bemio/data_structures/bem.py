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

'''
This class contains a structure to store hydrodynamic data from WAMTI,
AQWA, Nemoh, or another code that calculates hydrodynamic coefficients
and excitation forces

Author: Michael Lawson, Yi-Hsiang Yu, Carlos Michelen
'''
from __future__ import division

import numpy as np

import os

from scipy import interpolate

from scipy.linalg import hankel, expm

from progressbar import ProgressBar, Bar, Percentage

class Raw(object):
    '''
    Empty class to store raw data
    '''
    def __init__(self):
        pass

class HydrodynamicCoefficients(object):
    '''Hydrodynamic coefficients
    '''
    def __init__(self):
        self.irf             = ImpulseResponseFunction()
        self.ss              = StateSpaceRealization()

class ImpulseResponseFunction(object):
    '''Impulse response function data
    '''
    def __init__(self):
        pass

class StateSpaceRealization(object):
    '''State space realization data
    '''
    def __init__(self):
        pass

class HydrodynamicData(object):
    '''Hydrodynamic data from BEM simulations
    '''

    def __init__(self):
        # Default values
        self.rho            = 1000.         # Water density
        self.g              = 9.81          # Gravity
        self.wave_dir       = 0             # Wave directions
        self.num_bodies     = 0             # Number of bodies in the simulation
                     
        # np.array([])     
        self.cg             = 'not_defined' # Center of gravity                         
        self.cb             = 'not_defined' # Center of buoyancy                          
        self.k              = 'not_defined' # Linear hydrostatic spring stiffness                         
        self.T              = 'not_defined' # Wave periods                                                       
        self.w              = 'not_defined' # Wave frequencies                    
        
        # np.floats()
        self.wp_area        = 'not_defined'                             
        self.buoy_force     = 'not_defined'     
        self.disp_vol       = 'not_defined'                         
        self.water_depth    = 'not_defined'                           
        self.body_num       = 'not_defined'                    
        
        # strings
        self.dimensional        = 'not_defined' # Flag to determine if the hydrodynamic coefficients dimensional
        self.dimensionalize     = 'not_defined' # Flag to determine if the hydrodynamic coefficients should be output in dimensional or nondimensional form
        self.name               = 'not_defined' # Body name
        self.bem_code           = 'not_defined' # Name of the code used to produce the hydrodynamic coefficients
        self.bem_raw_data       = 'not_defined' # Raw data read by bemio

        # objects
        self.am              = HydrodynamicCoefficients()    
        self.rd              = HydrodynamicCoefficients()  
        self.ex              = HydrodynamicCoefficients()  
        self.rao             = HydrodynamicCoefficients()
        self.ssy             = HydrodynamicCoefficients()
        
    def __repr__(self):
        '''Custom output
        '''
        out_string = 'Body name: ' + str(self.name) + \
            '\n    Body number: ' + str(self.body_num) +\
            '\n    Total number of bodies: ' + str(self.num_bodies) + \
            '\n    Displaced volume (m^3): ' + str(self.disp_vol) + \
            '\n    Center of gravity (m): ' + str(self.cg) + \
            '\n    Center of buoyancy (m): ' + str(self.cb)
        return out_string


    def calc_irf_excitation(self, t_end=100.0, n_t = 1001, n_w=1001):
        '''Function to calculate the excitation impulse response function
        '''
        self.ex.irf.t = np.linspace(-t_end, t_end, n_t)
        self.ex.irf.w = np.linspace(self.w.min(),self.w.max(),n_w)

        self.ex.irf.f = np.zeros([self.ex.mag.shape[0], self.ex.mag.shape[1], self.ex.irf.t.size])

        ex_re_interp = interpolate_for_irf(self.w,self.ex.irf.w,self.ex.re)
        ex_im_interp = interpolate_for_irf(self.w,self.ex.irf.w,self.ex.im)

        pbar_maxval = self.ex.irf.t.size*self.ex.mag.shape[0]*self.ex.mag.shape[1]
        pbar = ProgressBar(widgets=['Calculating the excitation force impulse response function for ' + self.name + ':',Percentage(), Bar()], maxval=pbar_maxval).start()
        count = 1
        for t_ind, t in enumerate(self.ex.irf.t):

            for i in xrange(self.ex.mag.shape[0]):

                for j in xrange(self.ex.mag.shape[1]):
                    tmp = ex_re_interp[i,j,:]*np.cos(self.ex.irf.w*t) - ex_im_interp[i,j,:]*np.sin(self.ex.irf.w*t)
                    tmp *= 1./(2.*np.pi)
                    self.ex.irf.f[i,j,t_ind] = np.trapz(y=tmp,x=self.ex.irf.w)
                    pbar.update(count)
                    count += 1

        pbar.finish()

    def calc_ss_excitation(self):
        raise Exception('The calc_ss_excitation function is not yet implemented')

    def calc_irf_radiation(self, t_end=100, n_t = 1001, n_w=1001):
        '''Function to calculate the wave radiation impulse response function
        '''
        
        self.rd.irf.t = np.linspace(0,t_end,n_t)
        self.rd.irf.w = np.linspace(self.w.min(),self.w.max(),n_w)

        self.rd.irf.L = np.zeros( [ self.am.inf.shape[0],self.am.inf.shape[1],self.rd.irf.t.size ] )
        self.rd.irf.K = np.zeros( [ self.am.inf.shape[0],self.am.inf.shape[1],self.rd.irf.t.size ] )

        rd_interp = interpolate_for_irf(self.w,self.rd.irf.w,self.rd.all)

        # Calculate the IRF
        pbar = ProgressBar(widgets=['Calculating the radiation damping impulse response function for ' + self.name + ':',Percentage(), Bar()], maxval=np.size(self.rd.irf.t)*self.rd.all.shape[0]*self.rd.all.shape[1]).start()
        count = 1
        for t_ind, t in enumerate(self.rd.irf.t):

            for i in xrange(self.rd.all.shape[0]):

                for j in xrange(self.rd.all.shape[1]):
                    # Radiation damping calculation method
                    tmpL = 2./np.pi*rd_interp[i,j,:]*np.sin(self.rd.irf.w*t)
                    tmpK = 2./np.pi*rd_interp[i,j,:]*np.cos(self.rd.irf.w*t)

                    # Different IRF calculation methods are needed for dimensional and nondimensional hydro coefficients 
                    if self.dimensional is False:

                        tmpK *= self.rd.irf.w

                    elif self.dimensional is True:

                        tmpL /= self.rd.irf.w

                    self.rd.irf.K[i,j,t_ind] = np.trapz(y=tmpK,x=self.rd.irf.w)
                    self.rd.irf.L[i,j,t_ind] = np.trapz(y=tmpL,x=self.rd.irf.w)

                    pbar.update(count)
                    count += 1

        pbar.finish()

    def calc_ss_radiation(self, max_order=10, r2_thresh=0.95 ):
        '''Function to calculate state space realization
        
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
        dt                  = self.rd.irf.t[2]-self.rd.irf.t[1]
        numFreq             = np.size(self.rd.irf.t) 
        r2bt                = np.zeros( [ self.am.inf.shape[0],self.am.inf.shape[0],numFreq] )
        k_ss_est            = np.zeros( numFreq )
        self.rd.ss.irk_bss     = np.zeros( [ self.am.inf.shape[0],self.am.inf.shape[0],numFreq] )
        self.rd.ss.A           = np.zeros([6,self.am.inf.shape[1],max_order,max_order])
        self.rd.ss.B           = np.zeros([6,self.am.inf.shape[1],max_order,1])
        self.rd.ss.C           = np.zeros([6,self.am.inf.shape[1],1,max_order])
        self.rd.ss.D           = np.zeros([6,self.am.inf.shape[1],1])
        self.rd.ss.irk_bss     = np.zeros([6,self.am.inf.shape[1],numFreq])
        self.rd.ss.rad_conv    = np.zeros([6,self.am.inf.shape[1]])
        self.rd.ss.it          = np.zeros([6,self.am.inf.shape[1]])
        self.rd.ss.r2t         = np.zeros([6,self.am.inf.shape[1]])
        
        pbar = ProgressBar(widgets=['Calculating radiation damping state space coefficients for ' + self.name + ':',Percentage(), Bar()], maxval=self.am.inf.shape[0]*self.am.inf.shape[1]).start()
        count = 0
        for i in xrange(self.am.inf.shape[0]):

            for j in xrange(self.am.inf.shape[1]):

                r2bt = np.linalg.norm(self.rd.irf.K[i,j,:]-self.rd.irf.K.mean(axis=2)[i,j])
                
                ss = 2 #Initial state space order

                if r2bt != 0.0:
                    while True:
                        
                        #Perform Hankel Singular Value Decomposition
                        y=dt*self.rd.irf.K[i,j,:]                    
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
                        
                        ac = np.dot(CoeB*a-CoeD*np.eye(ss),iidd)                   #(A-I)2/T(I + A)^{-1}         = 2/T(A-I)(I + A)^{-1}
                        bc = (CoeA*CoeB-CoeC*CoeD)*np.dot(iidd,b)                  #(T/2+T/2)*2/T(I + A)^{-1}B   = 2(I + A)^{-1}B
                        cc = np.dot(c,iidd)                                        #C * 2/T(I + A)^{-1}          = 2/T(I + A)^{-1}
                        dc = d + CoeC*np.dot(np.dot(c,iidd),b)                     #D - T/2C (2/T(I + A)^{-1})B  = D - C(I + A)^{-1})B

                        for jj in xrange(numFreq):

                            k_ss_est[jj] = np.dot(np.dot(cc,expm(ac*dt*jj)),bc)    #Calculate impulse response function from state space approximation
      
                        R2TT = np.linalg.norm(self.rd.irf.K[i,j,:]-k_ss_est)          #Calculate 2 norm of the difference between know and estimated values impulse response function
                        R2T = 1 - np.square(R2TT/r2bt)                             #Calculate the R2 value for impulse response function

                        if R2T >= r2_thresh:                                       #Check to see if threshold for the impulse response is meet
                        
                            status = 1                                             #%Set status
                            break
                        
                        if ss == max_order:                                        #Check to see if limit on the state space order has been reached
                        
                            status = 2                                             #%Set status
                            break
                        
                        ss=ss+1                                                    #Increase state space order
                                            
                    self.rd.ss.A[i,j,0:ac.shape[0],0:ac.shape[0]]  = ac
                    self.rd.ss.B[i,j,0:bc.shape[0],0                ]  = bc[:,0]
                    self.rd.ss.C[i,j,0                ,0:cc.shape[1]]  = cc[0,:]
                    self.rd.ss.D[i,j]                                      = dc
                    self.rd.ss.irk_bss[i,j,:]  = k_ss_est
                    self.rd.ss.rad_conv[i,j] = status
                    self.rd.ss.r2t[i,j] = R2T
                    self.rd.ss.it[i,j] = ss

                count += 1
                pbar.update(count)

        pbar.finish()
        
  
    def dimensionalize_nondimensionalize(self):

        if self.dimensionalize is True and self.dimensional is False:

            self.k = self.k*self.rho*self.g
            self.am.all = self.am.all*self.rho
            self.am.inf = self.am.inf*self.rho
            self.am.zero = self.am.zero*self.rho
            self.ex.mag = self.ex.mag*self.rho*self.g
            self.ex.im = self.ex.im*self.rho*self.g
            self.ex.re = self.ex.re*self.rho*self.g

            for j in xrange(self.rd.all.shape[2]):

                self.rd.all[:,:,j] = self.rd.all[:,:,j]*self.rho*self.w[j]

            self.dimensional = True

            print 'Dimensionalizing hydro coefficients...'

            print 'Added mass, radiation damping, hydrodynamic excitation, and spring stiffness coefficients are dimensional'

        elif self.dimensionalize is False and self.dimensional is True:

            self.k = self.k/(self.rho*self.g) 
            self.am.all = self.am.all/self.rho
            self.am.inf = self.am.inf/self.rho
            self.am.zero = self.am.zero/self.rho
            self.ex.mag = self.ex.mag/(self.rho*self.g)
            self.ex.im = self.ex.im/(self.rho*self.g)
            self.ex.re = self.ex.re/(self.rho*self.g)

            for j in xrange(self.rd.all.shape[2]):

                    self.rd.all[:,:,j] = self.rd.all[:,:,j]/(self.rho*self.w[j])

            self.dimensional = False

            print 'Nondimesionailzing hydro coefficients...'

            print 'Added mass, radiation damping, hydrodynamic excitation, and spring stiffness coefficients are nondimensional'
            

def interpolate_for_irf(w_orig,w_interp,mat_in):
    '''
    Interpolate matrices for the IRF calculations
    '''
    mat_interp = np.zeros( [ mat_in.shape[0], mat_in.shape[1], w_interp.size ])

    flip = False

    if w_orig[0] > w_orig[1]:

        w_tmp = np.flipud(w_orig)
        flip = True

    else:

        w_tmp = w_orig

    for i in xrange(mat_in.shape[0]):

        for j in xrange(mat_in.shape[1]):

            if flip is True:

                rdTmp = np.flipud(mat_in[i,j,:])

            else:
                rdTmp = mat_in[i,j,:]

            f = interpolate.interp1d(x=w_tmp, y=rdTmp)
            mat_interp[i,j,:] = f(w_interp)

    return mat_interp

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

