# Copyright 2014 the National Renewable Energy Laboratory and Sandia Corporation

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import numpy as np

import os

from scipy import interpolate

from scipy.linalg import hankel, expm

from progressbar import ProgressBar, Bar, Percentage


class Raw(object):
    '''bemio.io.bem internal calss to store various data
    '''
    def __init__(self):
        pass


class HydrodynamicCoefficients(object):
    '''bemio.io.bem internal calss to store hydrodynamic coefficients
    '''
    def __init__(self):

        self.irf = ImpulseResponseFunction()
        self.ss = StateSpaceRealization()


class ImpulseResponseFunction(object):
    '''bemio.io.bem internal calss to store impulse response function data
    '''
    def __init__(self):
        pass


class StateSpaceRealization(object):
    '''bemio.io.bem internal calss to store state space realization data
    '''
    def __init__(self):
        pass


class HydrodynamicData(object):
    '''Class to store hydrodynamic data from NEMOH, WAMIT, and AQWA simulations

    Attribuites:
        rho : float
            Water density
        g : float
            Acceleration due to gravity
        wave_dir : np.array
            Array of wave directions used to calculate hydrodynamic coefficients
        num_bodies : int
            Number of bodies in the BEM simulation
        cg : np.array, shape = [3,1]
            Center of gravity
        cg : np.array, shape = [3,1]
            Center of buoyancy
        k : np.array, shape = [6,6]
            Linear hydrostatic restoring stiffness
        T : np.array
            Wave periods considered in BEM simulation
        w : np.array
            Wave frequencies considered in BEM simulation
        wp_area : float
            Water plan area
        buoy_force : float
            Buoyancy force acting on the body in the equelibrium position
        disp_volume : float
            Water volume displaced by the body in the equelibrium position
        water_depth : float
            Water depth considered in the BEM simulation
        body_num : int
            Body number of the body in the BEM simulation
        scaled : bool
            Flag that indicates if the hydrodynamic coefficients have been
            scaled. See the bemio.data_structures.bem.scale for more
            information.
        name : string
            Name of the body in the BEM simulation
        bem_code : string
            Name of the BWM code used to generate the hydrodynamic coefficients.
        bem_raw_data : string
            BEM output file
        am : HydrodynamicCoefficients
            Object containing the added mass coefficients. Attribuites of the am
            object are:
                am.all : np.array
                    Three dimensional array containing the added mass data.
                    The array is comprised of added mass matricies for each
                    frequency (`w`) of the BEM simulation. `am.all[0,0,:]` contains
                    the added mass matrix corresponding to `w[0]`, `am.all[1,1,:]`
                    contains the added mass matrix corresponding to `w[1]`, and so
                    on.
                am.inf : np.array
                    Two dimensional array containing the infiniate frequency added
                    mass data.
                am.zero : np.array
                    Two dimensional array containing the zero frequency added
                    mass data.
        ex : HydrodynamicCoefficients
            Object containing excitation force coefficinets. Attribuites of the
            ex object are:
                ex.mag : np.array
                    Magnitude of the excitation force for each frequency (`w`)
                    of the BEM simulation. `ex.re[0,:,:]` contains
                    the coefficients corresponding to `w[0]`, `ex.mag[1,:,:]`
                    contains the coefficients corresponding to `w[1]`, and so
                    on.
                ex.phase : np.array
                    Phase angle of the excitation force magnitude. The data
                    structure is organized the same way as ex.mag.
                ex.re : np.array
                    Real component of the excitation force. The data structure is
                    organized the same way as ex.mag.
                ex.im : np.array
                    Imaginary component of the excitation force. The data structure is
                    organized the same way as `ex.mag`.
                ex.sc : HydrodynamicCoefficients
                    A data object that contains wave scattering force coefficnets.
                    The orginization of the `ex.sc` coefficients is the same as the
                    `ex` object.
                ex.fk : HydrodynamicCoefficients
                    A data object that contains Froud-Krylof force coefficnets.
                    The orginization of the `ex.fk` coefficients is the same as the
                    `ex` object.
                ex.rao : HydrodynamicCoefficients
                    Object that contains response amplitude operator information.
                ex.ssy : HydrodynamicCoefficients
                    Object that contains WAMIT specific SSY data. See the WAMIT
                    users manual for more information. (Note: This variable
                    description needs to be improved)

    Note:
        This object is created by the bemio.io data readers and its purpose
        is to contain hydrodynamic coefficients from NEMOH, WAMIT, and AQWA
        data is a standard format.
    '''

    def __init__(self):
        # Default values
        self.rho = 1000.
        self.g = 9.81
        self.wave_dir = 0
        self.num_bodies = 0

        # np.array([])
        self.cg = 'not_defined'
        self.cb = 'not_defined'
        self.k = 'not_defined'
        self.T = 'not_defined'
        self.w = 'not_defined'

        # np.floats()
        self.wp_area = 'not_defined'
        self.buoy_force = 'not_defined'
        self.disp_vol = 'not_defined'
        self.water_depth = 'not_defined'
        self.body_num = 'not_defined'

        # strings`
        self.scaled = 'not_defined'
        self.name = 'not_defined'
        self.bem_code = 'not_defined'
        self.bem_raw_data = 'not_defined'

        # objects
        self.am = HydrodynamicCoefficients()
        self.rd = HydrodynamicCoefficients()
        self.ex = HydrodynamicCoefficients()
        self.ex.sc = HydrodynamicCoefficients()
        self.ex.fk = HydrodynamicCoefficients()
        self.rao = HydrodynamicCoefficients()
        self.ssy = HydrodynamicCoefficients()

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

    def calc_irf_excitation(self, t_end=100.0, n_t=1001, n_w=1001, w_min=None, w_max=None):
        '''Function to calculate the excitation impulse response function. For
        more information please see equation 2.12 in:

            Dynamics Modeling and Loads Analysis of an Offshore Floating Wind
            Turbine, J. M Jonkman, 2007, NREL/TP-500-41958

        Args:
            t_end : float
                Calculation range for the IRF. The IRF is calculated from -t_end
                to t_end
            n_t : int
                Number of timesteps in the IRF
            n_w : int
                Number of frequecy steps to use in the IRF calculation. If the
                frequency steps are different from the loaded BEM data, the
                excitation coefficinets are interpolated.

        Returns:
            No variables are directily returned by thi function. The following
            internal variables are calculated:

            self.ex.irf.t : np.array
                Timesteps at which the IRF is calculated
            self.ex.irf.w : np.array
                Frequency steps used to calculate the IRF
            self.ex.irf.f : np.array
                The excitation force IRF. This is a
            3 dimensional np.array:
                dim 1: number of degrees of freedom of the floating body,
                typically 6, dim 2: 1, dim 3: size(self.ex.irf.t)

        Note:
            When using this function it is important to look at the resulting
            IRF to insure that it is physical. The user may need to adjust the
            t_end, n_t, and n_w values in order to generate a realistic IRF.

        Examples:
            Once a HydrodynamicData data object is created (called 'hydro_data'
            here), the calc_irf_excitation function can be called:

            >>> hydro_data.calc_irf_excitation()

            The data can be plotted using matplotlib

            >>> import matplotlib.pyplot as plt
            >>> plt.plot(hydro_data.ex.irf.t,hydro_data.ex.irf.f)
            >>> plt.show()
        '''
        if w_min is None:
            w_min = self.w.min()
        if w_max is None:
            w_max = self.w.max()

        self.ex.irf.t = np.linspace(-t_end, t_end, n_t)
        self.ex.irf.w = np.linspace(w_min, w_max, n_w)

        self.ex.irf.f = np.zeros([self.ex.mag.shape[0], self.ex.mag.shape[1], self.ex.irf.t.size])

        ex_re_interp = _interpolate_for_irf(self.w, self.ex.irf.w, self.ex.re)
        ex_im_interp = _interpolate_for_irf(self.w, self.ex.irf.w, self.ex.im)

        pbar_maxval = self.ex.irf.t.size * self.ex.mag.shape[0] * self.ex.mag.shape[1]
        pbar = ProgressBar(widgets=['Excitation force IRF for ' + self.name + ':', Percentage(), Bar()], maxval=pbar_maxval).start()
        count = 1
        for t_ind, t in enumerate(self.ex.irf.t):

          for i in range(self.ex.mag.shape[0]):

            for j in range(self.ex.mag.shape[1]):
              tmp = ex_re_interp[i, j, :] * np.cos(self.ex.irf.w * t) - ex_im_interp[i, j, :] * np.sin(self.ex.irf.w * t)
              tmp *= 1. / np.pi
              self.ex.irf.f[i, j, t_ind] = np.trapz(y=tmp, x=self.ex.irf.w)
              pbar.update(count)
              count += 1

        pbar.finish()


    def calc_irf_radiation(self, t_end=100., n_t=1001, n_w=1001, w_min=None, w_max=None):
        '''Function to calculate the wave radiation impulse response function.
        For more information please see Section 2.4.2 in:

            Dynamics Modeling and Loads Analysis of an Offshore Floating Wind
            Turbine, J. M Jonkman, 2007, NREL/TP-500-41958,
            http://www.nrel.gov/docs/fy08osti/41958.pdf

        and Section 13.6-7 in:

            WAMIT v7.0 Users Manual, http://www.wamit.com/manual.htm


        Args:
            t_end : float
                Calculation range for the IRF. The IRF is calculated from 0 to
                t_end.
            n_t : int
                Number of timesteps in the IRF
            n_w : int
                Number of frequecy steps to use in the IRF calculation. If the
                frequency steps are different from the loaded BEM data, the
                radiation damping coefficinets are interpolated.
            w_min : float, optional
                Minimum frequency to use in IRF calculation. If no vlaue is
                given, this defaults to the minimum frequency from the BEM
                simulations
            w_max : float, optional
                Maximum frequency to use in IRF calculation. If no vlaue is
                given, this defaults to the maximum      frequency from the BEM
                simulations

        Returns:
            No variables are directily returned by thi function. The following
            internal variables are calculated:

            self.rd.irf.t : np.array
                Timesteps at which the IRF is calculated
            self.rd.irf.w : np.array
                Frequency steps used to calculate the IRF
            self.rd.irf.L : np.array
                The wave radiation IRF. This is a 3 dimensional np.array
                with dim 1: Number of degrees of freedom for the given body,
                dim 2: Number of degrees of freedom in the entire floating
                system. For example, two 6 degree of freedom bodies will
                result in a dimension of 12, dim 3: size(self.ex.irf.t)
            self.rd.irf.L : np.array
                d/dt(self.rd.irf.L)

        Note:
            When using this function it is important to look at the resulting
            IRF to insure that it is physical. The user may need to adjust the
            t_end, n_t, and n_w values in order to generate a realistic IRF.
        '''
        if w_min is None:
            w_min = self.w.min()
        if w_max is None:
            w_max = self.w.max()

        self.rd.irf.t = np.linspace(0, t_end, n_t)
        self.rd.irf.w = np.linspace(w_min, w_max, n_w)

        self.rd.irf.L = np.zeros(
            [self.am.inf.shape[0], self.am.inf.shape[1], self.rd.irf.t.size])
        self.rd.irf.K = np.zeros(
            [self.am.inf.shape[0], self.am.inf.shape[1], self.rd.irf.t.size])

        rd_interp = _interpolate_for_irf(self.w, self.rd.irf.w, self.rd.all)

        # Calculate the IRF
        pbar_max_val = self.rd.irf.t.size * self.rd.all.shape[0] * self.rd.all.shape[1]
        pbar = ProgressBar(widgets=['Radiation damping IRF for ' + self.name + ':', Percentage(), Bar()], maxval=pbar_max_val).start()
        count = 1
        for t_ind, t in enumerate(self.rd.irf.t):

            for i in range(self.rd.all.shape[0]):

                for j in range(self.rd.all.shape[1]):
              # Radiation damping calculation method
                    tmpL = 2. / np.pi * rd_interp[i, j, :] * np.sin(self.rd.irf.w * t)
                    tmpK = 2. / np.pi * rd_interp[i, j, :] * np.cos(self.rd.irf.w * t)

                    # Different IRF calculation methods are needed for dimensional and
                    # nondimensional hydro coefficients
                    if self.scaled is False:

                        tmpK *= self.rd.irf.w

                    elif self.scaled is True:

                        tmpL /= self.rd.irf.w

                    self.rd.irf.K[i, j, t_ind] = np.trapz(y=tmpK, x=self.rd.irf.w)
                    self.rd.irf.L[i, j, t_ind] = np.trapz(y=tmpL, x=self.rd.irf.w)

                    pbar.update(count)
                    count += 1

        pbar.finish()

    def calc_ss_radiation(self, max_order=10, r2_thresh=0.95):
        '''Function to calculate the state space reailization of the wave
        radiation IRF.

        Args:
            max_order : int
                Maximum order of the state space reailization fit
            r2_thresh : float
                The R^2 threshold used for the fit. THe value must be between 0
                and 1

        Returns:
            No variables are directily returned by thi function. The following
            internal variables are calculated:

            Ass :
                time-invariant state matrix
            Bss :
                time-invariant input matrix
            Css :
                time-invariant output matrix
            Dss :
                time-invariant feedthrough matrix
            k_ss_est :
                Impusle response function as cacluated from state space
                approximation
            status :
                status of the realization, 0 - zero hydrodynamic oefficients, 1
                - state space realization meets R2 thresholdm 2 - state space
                realization does not meet R2 threshold and at ss_max limit

        Examples:
        '''
        dt = self.rd.irf.t[2] - self.rd.irf.t[1]
        r2bt = np.zeros([self.am.inf.shape[0], self.am.inf.shape[0], self.rd.irf.t.size])
        k_ss_est = np.zeros(self.rd.irf.t.size)
        self.rd.ss.irk_bss = np.zeros([self.am.inf.shape[0], self.am.inf.shape[0], self.rd.irf.t.size])
        self.rd.ss.A = np.zeros([6, self.am.inf.shape[1], max_order, max_order])
        self.rd.ss.B = np.zeros([6, self.am.inf.shape[1], max_order, 1])
        self.rd.ss.C = np.zeros([6, self.am.inf.shape[1], 1, max_order])
        self.rd.ss.D = np.zeros([6, self.am.inf.shape[1], 1])
        self.rd.ss.irk_bss = np.zeros([6, self.am.inf.shape[1], self.rd.irf.t.size])
        self.rd.ss.rad_conv = np.zeros([6, self.am.inf.shape[1]])
        self.rd.ss.it = np.zeros([6, self.am.inf.shape[1]])
        self.rd.ss.r2t = np.zeros([6, self.am.inf.shape[1]])

        pbar = ProgressBar(widgets=['Radiation damping state space realization for ' + self.name + ':', Percentage(), Bar()], maxval=self.am.inf.shape[0] * self.am.inf.shape[1]).start()
        count = 0
        for i in range(self.am.inf.shape[0]):

          for j in range(self.am.inf.shape[1]):

            r2bt = np.linalg.norm(
                self.rd.irf.K[i, j, :] - self.rd.irf.K.mean(axis=2)[i, j])

            ss = 2  # Initial state space order

            if r2bt != 0.0:
              while True:

                # Perform Hankel Singular Value Decomposition
                y = dt * self.rd.irf.K[i, j, :]
                h = hankel(y[1::])
                u, svh, v = np.linalg.svd(h)

                u1 = u[0:self.rd.irf.t.size - 2, 0:ss]
                v1 = v.T[0:self.rd.irf.t.size - 2, 0:ss]
                u2 = u[1:self.rd.irf.t.size - 1, 0:ss]
                sqs = np.sqrt(svh[0:ss].reshape(ss, 1))
                invss = 1 / sqs
                ubar = np.dot(u1.T, u2)

                a = ubar * np.dot(invss, sqs.T)
                b = v1[0, :].reshape(ss, 1) * sqs
                c = u1[0, :].reshape(1, ss) * sqs.T
                d = y[0]

                CoeA = dt / 2
                CoeB = 1
                CoeC = -CoeA
                CoeD = 1

                # (T/2*I + T/2*A)^{-1}         = 2/T(I + A)^{-1}
                iidd = np.linalg.inv(CoeA * np.eye(ss) - CoeC * a)

                # (A-I)2/T(I + A)^{-1}         = 2/T(A-I)(I + A)^{-1}
                ac = np.dot(CoeB * a - CoeD * np.eye(ss), iidd)
                # (T/2+T/2)*2/T(I + A)^{-1}B   = 2(I + A)^{-1}B
                bc = (CoeA * CoeB - CoeC * CoeD) * np.dot(iidd, b)
                # C * 2/T(I + A)^{-1}          = 2/T(I + A)^{-1}
                cc = np.dot(c, iidd)
                # D - T/2C (2/T(I + A)^{-1})B  = D - C(I + A)^{-1})B
                dc = d + CoeC * np.dot(np.dot(c, iidd), b)

                for jj in range(self.rd.irf.t.size):

                  # Calculate impulse response function from state space
                  # approximation
                  k_ss_est[jj] = np.dot(np.dot(cc, expm(ac * dt * jj)), bc)

                # Calculate 2 norm of the difference between know and estimated
                # values impulse response function
                R2TT = np.linalg.norm(self.rd.irf.K[i, j, :] - k_ss_est)
                # Calculate the R2 value for impulse response function
                R2T = 1 - np.square(R2TT / r2bt)

                # Check to see if threshold for the impulse response is meet
                if R2T >= r2_thresh:

                  status = 1  # %Set status
                  break

                # Check to see if limit on the state space order has been reached
                if ss == max_order:

                  status = 2  # %Set status
                  break

                ss = ss + 1  # Increase state space order

              self.rd.ss.A[i, j, 0:ac.shape[0], 0:ac.shape[0]] = ac
              self.rd.ss.B[i, j, 0:bc.shape[0], 0] = bc[:, 0]
              self.rd.ss.C[i, j, 0, 0:cc.shape[1]] = cc[0, :]
              self.rd.ss.D[i, j] = dc
              self.rd.ss.irk_bss[i, j, :] = k_ss_est
              self.rd.ss.rad_conv[i, j] = status
              self.rd.ss.r2t[i, j] = R2T
              self.rd.ss.it[i, j] = ss

            count += 1
            pbar.update(count)

        pbar.finish()

    def scale(self, scale=None):
        '''Function to scale the hydrodynamic coefficient.

        Args:
            scale (bool): Boolean operater to determin if hydrodynamic data
                should be scaled. If `scale is True` self.am (added mass)
                coeffcients are scaled by self.rho*self.g, self.ex (excitation)
                coefficinets are scaled by self.rho*self.g, and self.rd
                (radiation damping) coefficnets are scaled by self.rho*self.w.
                If `scale is False` self.am (added mass) coeffcients are scaled by
                1./(self.rho*self.g), self.ex (excitation) coefficinets are scaled by
                1./(self.rho*self.g), and self.rd (radiation damping) coefficnets are
                scaled by 1./(self.rho*self.w).

        Returns:
            No variables are directily returned by thi function. Hydrodynamic
            coefficnets are scaled as described above.

        Note:
            The bemio.io.nemoh, bemio.io.wamit, and bemio.io.aqwa functions read
            data and return it with the `scale == False`. If the user wishes to
            return data from these functions with `scale == True`, they shoud
            specify the `scale` argument when the data is read or should call
            the scale function after read. Regardless, for consistency, the
            data should be scaled before the `calc_irf_radiation`,
            `calc_irf_excitation`, or `calc_ss_radiation` functions are called.

        Examples:
            Once a HydrodynamicData data object is created (called 'hydro_data'
            here), the scalen function can be called:

            >>> hydro_data.scale(scale=False)
        '''
        if scale is not None:
            self.scale = scale

        if self.scale is True and self.scaled is False:
          print('\tScaling hydro coefficients for body ' + self.name + ' by rho, g, and w...')
          try:
              self.k *= self.rho * self.g
          except:
              print('\t\tSpring stiffness not scaled')

          self.am.all *= self.rho
          self.am.inf *= self.rho
          if hasattr(self.am,'zero') is True:
              self.am.zero *= self.rho

          self.ex.mag *= self.rho * self.g
          self.ex.re *= self.rho * self.g
          self.ex.im *= self.rho * self.g

          if hasattr(self.ex.sc,'mag') is True:
              self.ex.sc.mag *= self.rho * self.g
              self.ex.sc.re *= self.rho * self.g
              self.ex.sc.im *= self.rho * self.g

          if hasattr(self.ex.fk,'mag') is True:
              self.ex.fk.mag *= self.rho * self.g
              self.ex.fk.re *= self.rho * self.g
              self.ex.fk.im *= self.rho * self.g


          for j in range(self.rd.all.shape[2]):

            self.rd.all[:, :, j] = self.rd.all[:, :, j] * self.rho * self.w[j]

          self.scaled = True

        elif self.scale is False and self.scaled is True:
          print('\tUn-scaling hydro coefficients for body ' + self.name + ' by rho, g, and w...')
          try:
              self.k /= (self.rho * self.g)
          except:
              print('\t\tSpring stiffness not un-scaled')
          self.am.all /= self.rho
          self.am.inf /= self.rho
          if hasattr(self.am,'zero') is True:
              self.am.zero /= self.rho

          self.ex.mag /= (self.rho * self.g)
          self.ex.re /= (self.rho * self.g)
          self.ex.im /= (self.rho * self.g)

          if hasattr(self.ex.sc,'mag') is True:
              self.ex.sc.mag /= (self.rho * self.g)
              self.ex.sc.re /= (self.rho * self.g)
              self.ex.sc.im /= (self.rho * self.g)

          if hasattr(self.ex.fk,'mag') is True:
              self.ex.fk.mag /= (self.rho * self.g)
              self.ex.fk.re /= (self.rho * self.g)
              self.ex.fk.im /= (self.rho * self.g)

          for j in range(self.rd.all.shape[2]):

            self.rd.all[:, :, j] = self.rd.all[:, :, j] / (self.rho * self.w[j])

          self.scaled = False


def _interpolate_for_irf(w_orig, w_interp, mat_in):
  '''
  Interpolate matrices for the IRF calculations
  '''
  mat_interp = np.zeros([mat_in.shape[0], mat_in.shape[1], w_interp.size])

  flip = False

  if w_orig[0] > w_orig[1]:

    w_tmp = np.flipud(w_orig)
    flip = True

  else:

    w_tmp = w_orig

  for i in range(mat_in.shape[0]):

    for j in range(mat_in.shape[1]):

      if flip is True:

        rdTmp = np.flipud(mat_in[i, j, :])

      else:
        rdTmp = mat_in[i, j, :]

      f = interpolate.interp1d(x=w_tmp, y=rdTmp)
      mat_interp[i, j, :] = f(w_interp)

  return mat_interp


def generate_file_names(out_file):
  '''
  Internal bemio function to generate filenames needed by hydroData module

  Parameters:
    out_file : str
        Name of hydrodynamic data file

  Returns:
    files : dictionary
        A dictionary of file names used by bemio
  '''
  out_file = os.path.abspath(out_file)
  (path, file_name) = os.path.split(out_file)
  file_name_raw, file_extension =  os.path.splitext(file_name)

  files = {}
  files['base_name'] = os.path.join(path, file_name_raw)
  files['out'] = os.path.join(path, file_name)
  files['hdf5'] = os.path.join(path, file_name_raw + '.h5')
  files['pickle'] = os.path.join(path, file_name_raw + '.p')

  return files
