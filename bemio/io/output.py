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
import h5py
import numpy as np
from bemio.__version__ import base

def write_hdf5(bemio_obj,out_file=None):
    '''
    Function that writes NEMOH, WAMIT, or NEMOH data to a standard human
    readable data format that uses the HDF5 format. This data can easily be
    input into various codes, such as MATLAB, Python, C++, etc. The data can
    easily be viewed using `HDFVIEW <https://www.hdfgroup.org/products/java/hdfview/>`_.

    Parameters:
        data_object : {bemio.io.NemohOutput, bamio.io.WamitOutput, bemio.io.AqwaOutput}
            A data object created using the bemio.io data readers
        out_file : str, optional
            The name of the output file. The file should have the .h5 file
            extension

    Examples:
        This example assumes there is a wamit output file called wamit.out that
        is read using the bemio.io.wamit.read function

        >>> from bemio.io.wamit import read
        >>> from bemio.io.output import write_hdf5
        >>> wamit_data = read(out_file=wamit.out)
        >>> write_hdf5(wamit_data)
        Writing HDF5 data to ./wamit.h5
    '''


    if out_file is None:
        out_file = bemio_obj.files['hdf5']


    print('Writing HDF5 data to ' + out_file)


    with h5py.File(out_file, "w") as f:

        for key, key in enumerate(bemio_obj.body.keys()):

            # Body properities
            cg = f.create_dataset('body' + str(key+1) + '/properties/cg',data=bemio_obj.body[key].cg)
            cg.attrs['units'] = 'm'
            cg.attrs['description'] = 'Center of gravity'

            cb = f.create_dataset('body' + str(key+1) + '/properties/cb',data=bemio_obj.body[key].cb)
            cb.attrs['units'] = 'm'
            cb.attrs['description'] = 'Center of buoyancy'

            vol = f.create_dataset('body' + str(key+1) + '/properties/disp_vol',data=bemio_obj.body[key].disp_vol)
            vol.attrs['units'] = 'm^3'
            vol.attrs['description'] = 'Displaced volume'

            name = f.create_dataset('body' + str(key+1) + '/properties/name',data=bemio_obj.body[key].name)
            name.attrs['description'] = 'Name of rigid body'

            num = f.create_dataset('body' + str(key+1) + '/properties/body_number',data=bemio_obj.body[key].body_num)
            num.attrs['description'] = 'Number of rigid body from the BEM simulation'

            # Hydro coeffs
            # Radiation IRF
            try:

                irf_rad_k_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/K',data=bemio_obj.body[key].rd.irf.K)
                irf_rad_k_correct_loc.attrs['units'] = ''
                irf_rad_k_correct_loc.attrs['description'] = 'Impulse response function'

                irf_rad_t_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/t',data=bemio_obj.body[key].rd.irf.t)
                irf_rad_t_correct_loc.attrs['units'] = 'seconds'
                irf_rad_t_correct_loc.attrs['description'] = 'Time vector for the impulse response function'

                irf_rad_w_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/w',data=bemio_obj.body[key].rd.irf.w)
                irf_rad_w_correct_loc.attrs['units'] = 'seconds'
                irf_rad_w_correct_loc.attrs['description'] = 'Interpolated frequencies used to compute the impulse response function'

                irf_rad_l_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/L',data=bemio_obj.body[key].rd.irf.L)
                irf_rad_l_correct_loc.attrs['units'] = ''
                irf_rad_l_correct_loc.attrs['description'] = 'Time derivative of the impulse response function'


                for m in range(bemio_obj.body[key].am.all.shape[0]):

                    for n in range(bemio_obj.body[key].am.all.shape[1]):

                        irf_rad_l_comp_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/components/L/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].rd.irf.t,bemio_obj.body[key].rd.irf.L[m,n,:]]).transpose())
                        irf_rad_l_comp_correct_loc.attrs['units'] = ''
                        irf_rad_l_comp_correct_loc.attrs['description'] = 'Components of the IRF'

                        irf_rad_k_comp_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/components/K/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].rd.irf.t,bemio_obj.body[key].rd.irf.K[m,n,:]]).transpose())
                        irf_rad_k_comp_correct_loc.attrs['units'] = ''
                        irf_rad_k_comp_correct_loc.attrs['description'] = 'Components of the ddt(IRF): K'
            except:

                print('\tRadiation IRF functions for ' + bemio_obj.body[key].name + ' were not written.')

            # Excitation IRF
            try:

                irf_ex_f = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/impulse_response_fun/f',data=bemio_obj.body[key].ex.irf.f)
                irf_ex_f.attrs['units'] = ''
                irf_ex_f.attrs['description'] = 'Impulse response function'

                irf_ex_t = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/impulse_response_fun/w',data=bemio_obj.body[key].ex.irf.w)
                irf_ex_w = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/impulse_response_fun/t',data=bemio_obj.body[key].ex.irf.t)

                for m in range(bemio_obj.body[key].ex.mag.shape[0]):

                    for n in range(bemio_obj.body[key].ex.mag.shape[1]):

                        irf_ex_f_comp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/impulse_response_fun/components/f/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].ex.irf.t,bemio_obj.body[key].ex.irf.f[m,n,:]]).transpose())
                        irf_ex_f_comp.attrs['units'] = ''
                        irf_ex_f_comp.attrs['description'] = 'Components of the ddt(IRF): f'

            except:

                print('\tExcitation IRF functions for ' + bemio_obj.body[key].name + ' were not written.')

            try:

                ssRadfA = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/A/all',data=bemio_obj.body[key].rd.ss.A)
                ssRadfA.attrs['units'] = ''
                ssRadfA.attrs['description'] = 'State Space A Coefficient'

                ssRadfB = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/B/all',data=bemio_obj.body[key].rd.ss.B)
                ssRadfB.attrs['units'] = ''
                ssRadfB.attrs['description'] = 'State Space B Coefficient'

                ssRadfC = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/C/all',data=bemio_obj.body[key].rd.ss.C)
                ssRadfC.attrs['units'] = ''
                ssRadfC.attrs['description'] = 'State Space C Coefficient'

                ssRadfD = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/D/all',data=bemio_obj.body[key].rd.ss.D)
                ssRadfD.attrs['units'] = ''
                ssRadfD.attrs['description'] = 'State Space D Coefficient'

                r2t = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/r2t',data=bemio_obj.body[key].rd.ss.r2t)
                r2t.attrs['units'] = ''
                r2t.attrs['description'] = 'State space curve fitting R**2 value'

                it = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/it',data=bemio_obj.body[key].rd.ss.it)
                it.attrs['units'] = ''
                it.attrs['description'] = 'Order of state space realization'

                for m in range(bemio_obj.body[key].am.all.shape[0]):

                    for n in range(bemio_obj.body[key].am.all.shape[1]):

                        ss_A = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/A/components/' + str(m+1) + '_' + str(n+1),data=bemio_obj.body[key].rd.ss.A[m,n,:,:])
                        ss_A.attrs['units'] = ''
                        ss_A.attrs['description'] = 'Components of the State Space A Coefficient'

                        ss_B = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/B/components/' + str(m+1) + '_' + str(n+1),data=bemio_obj.body[key].rd.ss.B[m,n,:,:])
                        ss_B.attrs['units'] = ''
                        ss_B.attrs['description'] = 'Components of the State Space B Coefficient'

                        ss_C = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/C/components/' + str(m+1) + '_' + str(n+1),data=bemio_obj.body[key].rd.ss.C[m,n,:,:])
                        ss_C.attrs['units'] = ''
                        ss_C.attrs['description'] = 'Components of the State Space C Coefficient'

                        ss_D = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/state_space/D/components/' + str(m+1) + '_' + str(n+1),data=bemio_obj.body[key].rd.ss.D[m,n])
                        ss_D.attrs['units'] = ''
                        ss_D.attrs['description'] = 'Components of the State Space C Coefficient'

            except:

                print('\tRadiation state space coefficients for ' + bemio_obj.body[key].name + ' were not written.')

            k = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/linear_restoring_stiffness',data=bemio_obj.body[key].k)
            k.attrs['units'] = ''
            k.attrs['description'] = 'Hydrostatic stiffness matrix'

            exMag = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/mag',data=bemio_obj.body[key].ex.mag)
            exMag.attrs['units'] = ''
            exMag.attrs['description'] = 'Magnitude of excitation force'

            exPhase = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/phase',data=bemio_obj.body[key].ex.phase)
            exPhase.attrs['units'] = 'rad'
            exPhase.attrs['description'] = 'Phase angle of excitation force'

            exRe = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/re',data=bemio_obj.body[key].ex.re)
            exRe.attrs['units'] = ''
            exRe.attrs['description'] = 'Real component of excitation force'

            exIm = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/im',data=bemio_obj.body[key].ex.im)
            exIm.attrs['units'] = ''
            exIm.attrs['description'] = 'Imaginary component of excitation force'

            for m in range(bemio_obj.body[key].ex.mag.shape[0]):

                for n in range(bemio_obj.body[key].ex.mag.shape[1]):

                    irf_ex_f_comp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation//components/mag/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].T,bemio_obj.body[key].ex.mag[m,n,:]]).transpose())
                    irf_ex_f_comp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation//components/phase/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].T,bemio_obj.body[key].ex.phase[m,n,:]]).transpose())
                    irf_ex_f_comp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation//components/re/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].T,bemio_obj.body[key].ex.re[m,n,:]]).transpose())
                    irf_ex_f_comp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation//components/im/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].T,bemio_obj.body[key].ex.im[m,n,:]]).transpose())


            # Scattering and FK forces
            try:
                ex_sc_Mag = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/scattering/mag',data=bemio_obj.body[key].ex.sc.mag)
                ex_sc_Mag.attrs['units'] = ''
                ex_sc_Mag.attrs['description'] = 'Magnitude of excitation force'

                ex_sc_Phase = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/scattering/phase',data=bemio_obj.body[key].ex.sc.phase)
                ex_sc_Phase.attrs['units'] = 'rad'
                ex_sc_Phase.attrs['description'] = 'Phase angle of excitation force'

                ex_sc_Re = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/scattering/re',data=bemio_obj.body[key].ex.sc.re)
                ex_sc_Re.attrs['units'] = ''
                ex_sc_Re.attrs['description'] = 'Real component of excitation force'

                ex_sc_Im = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/scattering/im',data=bemio_obj.body[key].ex.sc.im)
                ex_sc_Im.attrs['units'] = ''
                ex_sc_Im.attrs['description'] = 'Imaginary component of excitation force'

                ex_fk_Mag = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/froud_krylof/mag',data=bemio_obj.body[key].ex.fk.mag)
                ex_fk_Mag.attrs['units'] = ''
                ex_fk_Mag.attrs['description'] = 'Magnitude of excitation force'

                ex_fk_Phase = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/froud_krylof/phase',data=bemio_obj.body[key].ex.fk.phase)
                ex_fk_Phase.attrs['units'] = 'rad'
                ex_fk_Phase.attrs['description'] = 'Phase angle of excitation force'

                ex_fk_Re = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/froud_krylof/re',data=bemio_obj.body[key].ex.fk.re)
                ex_fk_Re.attrs['units'] = ''
                ex_fk_Re.attrs['description'] = 'Real component of excitation force'

                ex_fk_Im = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/froud_krylof/im',data=bemio_obj.body[key].ex.fk.im)
                ex_fk_Im.attrs['units'] = ''
                ex_fk_Im.attrs['description'] = 'Imaginary component of excitation force'

            except:
                pass

            # Write added mass information
            amInf = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/inf_freq',data=bemio_obj.body[key].am.inf)
            amInf.attrs['units for translational degrees of freedom'] = 'kg'
            amInf.attrs['description'] = 'Infinite frequency added mass'

            am = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/all',data=bemio_obj.body[key].am.all)
            am.attrs['units for translational degrees of freedom'] = 'kg'
            am.attrs['units for rotational degrees of freedom'] = 'kg-m^2'
            am.attrs['description'] = 'Added mass. Frequency is the third dimension of the data structure.'

            for m in range(bemio_obj.body[key].am.all.shape[0]):

                for n in range(bemio_obj.body[key].am.all.shape[1]):

                    amComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/components/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].T, bemio_obj.body[key].am.all[m,n,:]]).transpose())
                    amComp.attrs['units'] = ''
                    amComp.attrs['description'] = 'Added mass components as a function of frequency'

                    radComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/components/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].T, bemio_obj.body[key].rd.all[m,n,:]]).transpose())
                    radComp.attrs['units'] = ''
                    radComp.attrs['description'] = 'Radiation damping components as a function of frequency'

            rad = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/all',data=bemio_obj.body[key].rd.all)
            rad.attrs['units'] = ''
            rad.attrs['description'] = 'Radiation damping. Frequency is the thrid dimension of the data structure.'

        # Simulation parameters
        g = f.create_dataset('simulation_parameters/g',data=bemio_obj.body[key].g)
        g.attrs['units'] = 'm/s^2'
        g.attrs['description'] = 'Gravitational acceleration'

        rho = f.create_dataset('simulation_parameters/rho',data=bemio_obj.body[key].rho)
        rho.attrs['units'] = 'kg/m^3'
        rho.attrs['description'] = 'Water density'

        T = f.create_dataset('simulation_parameters/T',data=bemio_obj.body[key].T)
        T.attrs['units'] = 's'
        T.attrs['description'] = 'Wave periods'

        w = f.create_dataset('simulation_parameters/w',data=bemio_obj.body[key].w)
        w.attrs['units'] = 'rad/s'
        w.attrs['description'] = 'Wave frequencies'

        water_depth = f.create_dataset('simulation_parameters/water_depth',data=bemio_obj.body[key].water_depth)
        water_depth.attrs['units'] = 'm'
        water_depth.attrs['description'] = 'Water depth'

        wave_dir = f.create_dataset('simulation_parameters/wave_dir',data=bemio_obj.body[key].wave_dir)
        wave_dir.attrs['units'] = 'rad'
        wave_dir.attrs['description'] = 'Wave direction'

        scaled = f.create_dataset('simulation_parameters/scaled',data=bemio_obj.body[key].scaled)
        scaled.attrs['description'] = 'True: The data is scaled by rho*g, False: The data is not scaled by rho*g'

        bemio_version = f.create_dataset('bemio_information/version',data=base())

        rawOut = f.create_dataset('bem_data/output',data=bemio_obj.body[key].bem_raw_data)
        rawOut.attrs['description'] = 'Raw output from BEM code'

        code = f.create_dataset('bem_data/code',data=bemio_obj.body[key].bem_code)
        code.attrs['description'] = 'BEM code'
