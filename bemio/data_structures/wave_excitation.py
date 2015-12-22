# Copyright 2014 the National Renewable Energy Laboratory and Sandia Corporation
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import numpy as np
from scipy import interpolate
from progressbar import ProgressBar, Bar, Percentage

class ImpulseResponseFunction(object):
    '''Internal bemio object to contain impulse response function (IRF) data
    '''
    pass

class WaveElevationTimeSeries(object):
    '''Internal bemio object to contain wave elevation time series data
    '''
    pass

class WaveExcitationForce(object):
    '''Internal bemio object to contain wave excitation force data
    '''
    pass

class WaveExcitationConvolution(object):
    '''
    Object for calculating wave excitation force time history using the
    convolution method

    Parameters:
        irf : np.array
            Wave excitation force impulse response function.
        irf_t : np.array
            Time series corresponding to `irf`
        eta : np.array
            Wave elevation time series
        eta_t : np.array
            Time series corresponding to `eta`

    Attribuites:
        self.irf : ImpulseResponseFunction
            Object containing excitation force IRF information
        self.wave_elevation : WaveElevationTimeSeries
            Object containing wave elevation time series data
        self.excitation_force : WaveExcitationForce
            Object containing wave excitation force data
    '''
    def __init__(self, irf, irf_t, eta, eta_t):

        self.irf = ImpulseResponseFunction()
        self.wave_elevation = WaveElevationTimeSeries()
        self.excitation_force = WaveExcitationForce()

        self.irf.f = irf
        self.irf.t = irf_t

        self.wave_elevation.eta = eta
        self.wave_elevation.t = eta_t
        self.wave_elevation.dt = self.wave_elevation.t[1] - self.wave_elevation.t[0]
        self._excitation_convolution()

    def _excitation_convolution(self):
        '''Internal function to perform the wave excitation convolution
        '''
        eta_interp = interpolate.interp1d(x=self.wave_elevation.t, y=self.wave_elevation.eta, bounds_error=False, fill_value=0.)
        irf_interp = interpolate.interp1d(x=self.irf.t, y=self.irf.f, bounds_error=False, fill_value=0.)

        # Interpolate the IRF to the dt as the wave elevation data
        irf = irf_interp(np.linspace(self.irf.t.min(),self.irf.t.max(),(self.irf.t.max()-self.irf.t.min())/self.wave_elevation.dt+1))

        # Assume that the IRF dt is used unless specified by the user
        # if self.excitation_force.dt is None:
        #     self.excitation_force.dt = self.irf.t[1] - self.irf.t[0]
        # This code caluclates the wave excitation force manually - the below method that uses the convolve function is much more efficient
        # self.excitation_force.t = np.linspace(self.wave_elevation.t.min(), self.wave_elevation.t.max(), (self.wave_elevation.t.max()-self.wave_elevation.t.min())/self.excitation_force.dt+1)
        # pbar_max_val = self.excitation_force.t.max()
        # pbar = ProgressBar(widgets=['Calculating the excitation force time history:', Percentage(), Bar()], maxval=pbar_max_val).start()
        # f_ex = []
        # for t in self.excitation_force.t:
        #     f_ex.append(np.trapz(y=irf_interp(self.irf.t)*eta_interp(t-self.irf.t),x=self.irf.t))
        #
        #     pbar.update(t)
        # pbar.finish()

        f_ex_conv = np.convolve(self.wave_elevation.eta, irf, mode='same')*self.wave_elevation.dt

        self.excitation_force.f = np.array(f_ex_conv)
        self.excitation_force.t = self.wave_elevation.t


def convolution(irf, irf_t, eta, eta_t, dt=None):
    '''
    Function to calculate wave excitation force using the convolution method

    Patrameters:
        irf : np.array
            Wave excitation force impulse response function.
        irf_t : np.array
            Time series corresponding to `irf`
        eta : np.array
            Wave elevation time series
        eta_t : np.array
            Time series corresponding to `eta`
        dt : float, optional
            Time step for calculating the

    Returns:
        excitation_force : WaveExcitationConvolution
            This function returns a `WaveExcitationConvolution` object with
            the wave exciting force and other information. See the
            `WaveExcitationConvolution` for more information.

    Example:
        The following example assumes that variables `irf`, `irf_t`, `eta`, and
        `eta_t` of type type(np.array) exist in the workspace. The contents of
        these variables are described above.

        Calculate excitation force using the convolution method

        >>> ex = convolution(irf=irf, irf_t=irf_t, eta=eta, eta_t=eta_t)

        Plot the data

        >>> plt.figure()
        >>> plt.plot(ex.excitation_force.t,ex.excitation_force.f)
    '''
    excitation_force = WaveExcitationConvolution(irf, irf_t, eta, eta_t)
    return excitation_force
