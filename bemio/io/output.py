import h5py
import numpy as np
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


    print 'Writing HDF5 data to ' + out_file


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
            # Radiation IRF
            try:

                irf_rad_k = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/K',data=data_obj.data[key].rd.irf.K)
                irf_rad_k.attrs['units'] = ''
                irf_rad_k.attrs['description'] = 'Impulse response function'

                irf_rad_t = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/t',data=data_obj.data[key].rd.irf.t)
                irf_rad_t.attrs['units'] = 'seconds'
                irf_rad_t.attrs['description'] = 'Time vector for the impulse response function'

                irf_rad_w = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/w',data=data_obj.data[key].rd.irf.w)
                irf_rad_w.attrs['units'] = 'seconds'
                irf_rad_w.attrs['description'] = 'Interpolated frequencies used to compute the impulse response function'

                irf_rad_l = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/L',data=data_obj.data[key].rd.irf.L)
                irf_rad_l.attrs['units'] = ''
                irf_rad_l.attrs['description'] = 'Time derivative of the impulse response function'

                irf_rad_k_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/K',data=data_obj.data[key].rd.irf.K)
                irf_rad_k_correct_loc.attrs['units'] = ''
                irf_rad_k_correct_loc.attrs['description'] = 'Impulse response function'

                irf_rad_t_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/t',data=data_obj.data[key].rd.irf.t)
                irf_rad_t_correct_loc.attrs['units'] = 'seconds'
                irf_rad_t_correct_loc.attrs['description'] = 'Time vector for the impulse response function'

                irf_rad_w_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/w',data=data_obj.data[key].rd.irf.w)
                irf_rad_w_correct_loc.attrs['units'] = 'seconds'
                irf_rad_w_correct_loc.attrs['description'] = 'Interpolated frequencies used to compute the impulse response function'

                irf_rad_l_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/L',data=data_obj.data[key].rd.irf.L)
                irf_rad_l_correct_loc.attrs['units'] = ''
                irf_rad_l_correct_loc.attrs['description'] = 'Time derivative of the impulse response function'


                for m in xrange(data_obj.data[key].am.all.shape[0]):

                    for n in xrange(data_obj.data[key].am.all.shape[1]):

                        irf_rad_l_comp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/components/L/' + str(m+1) + '_' + str(n+1),data=np.array([data_obj.data[key].rd.irf.t,data_obj.data[key].rd.irf.L[m,n,:]]).transpose())
                        irf_rad_l_comp.attrs['units'] = ''
                        irf_rad_l_comp.attrs['description'] = 'Components of the IRF'

                        irf_rad_k_comp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/impulse_response_fun/components/K/' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].rd.irf.K[m,n,:])
                        irf_rad_k_comp.attrs['units'] = ''
                        irf_rad_k_comp.attrs['description'] = 'Components of the ddt(IRF): K'

                        irf_rad_l_comp_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/components/L/' + str(m+1) + '_' + str(n+1),data=np.array([data_obj.data[key].rd.irf.t,data_obj.data[key].rd.irf.L[m,n,:]]).transpose())
                        irf_rad_l_comp_correct_loc.attrs['units'] = ''
                        irf_rad_l_comp_correct_loc.attrs['description'] = 'Components of the IRF'

                        irf_rad_k_comp_correct_loc = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/impulse_response_fun/components/K/' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].rd.irf.K[m,n,:])
                        irf_rad_k_comp_correct_loc.attrs['units'] = ''
                        irf_rad_k_comp_correct_loc.attrs['description'] = 'Components of the ddt(IRF): K'
            except:

                print '\tRadiation IRF functions for ' + data_obj.data[key].name + ' were not written.'

            # Excitation IRF
            try:

                irf_ex_f = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/impulse_response_fun/f',data=data_obj.data[key].ex.irf.f)
                irf_ex_f.attrs['units'] = ''
                irf_ex_f.attrs['description'] = 'Impulse response function'

                irf_ex_t = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/impulse_response_fun/w',data=data_obj.data[key].ex.irf.w)
                irf_ex_w = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/impulse_response_fun/t',data=data_obj.data[key].ex.irf.t)

                for m in xrange(data_obj.data[key].ex.mag.shape[0]):

                    for n in xrange(data_obj.data[key].ex.mag.shape[1]):

                        irf_ex_f_comp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/impulse_response_fun/components/f/' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].ex.irf.f[m,n,:])
                        irf_ex_f_comp.attrs['units'] = ''
                        irf_ex_f_comp.attrs['description'] = 'Components of the ddt(IRF): f'

            except:

                print '\tExcitation IRF functions for ' + data_obj.data[key].name + ' were not written.'

            try:

                ssRadfA = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/A/all',data=data_obj.data[key].rd.ss.A)
                ssRadfA.attrs['units'] = ''
                ssRadfA.attrs['description'] = 'State Space A Coefficient'

                ssRadfB = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/B/all',data=data_obj.data[key].rd.ss.B)
                ssRadfB.attrs['units'] = ''
                ssRadfB.attrs['description'] = 'State Space B Coefficient'

                ssRadfC = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/C/all',data=data_obj.data[key].rd.ss.C)
                ssRadfC.attrs['units'] = ''
                ssRadfC.attrs['description'] = 'State Space C Coefficient'

                ssRadfD = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/D/all',data=data_obj.data[key].rd.ss.D)
                ssRadfD.attrs['units'] = ''
                ssRadfD.attrs['description'] = 'State Space D Coefficient'

                r2t = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/r2t',data=data_obj.data[key].rd.ss.r2t)
                r2t.attrs['units'] = ''
                r2t.attrs['description'] = 'State space curve fitting R**2 value'

                it = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/it',data=data_obj.data[key].rd.ss.it)
                it.attrs['units'] = ''
                it.attrs['description'] = 'Order of state space realization'

                for m in xrange(data_obj.data[key].am.all.shape[0]):

                    for n in xrange(data_obj.data[key].am.all.shape[1]):

                        ss_A = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/A/components/' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].rd.ss.A[m,n,:,:])
                        ss_A.attrs['units'] = ''
                        ss_A.attrs['description'] = 'Components of the State Space A Coefficient'

                        ss_B = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/B/components/' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].rd.ss.B[m,n,:,:])
                        ss_B.attrs['units'] = ''
                        ss_B.attrs['description'] = 'Components of the State Space B Coefficient'

                        ss_C = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/C/components/' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].rd.ss.C[m,n,:,:])
                        ss_C.attrs['units'] = ''
                        ss_C.attrs['description'] = 'Components of the State Space C Coefficient'

                        ss_D = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/state_space/D/components/' + str(m+1) + '_' + str(n+1),data=data_obj.data[key].rd.ss.D[m,n])
                        ss_D.attrs['units'] = ''
                        ss_D.attrs['description'] = 'Components of the State Space C Coefficient'

            except:

                print '\tRadiation state space coefficients for ' + data_obj.data[key].name + ' were not written.'

            k = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/linear_restoring_stiffness',data=data_obj.data[key].k)
            k.attrs['units'] = ''
            k.attrs['description'] = 'Hydrostatic stiffness matrix'

            exMag = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/mag',data=data_obj.data[key].ex.mag)
            exMag.attrs['units'] = ''
            exMag.attrs['description'] = 'Magnitude of excitation force'

            exPhase = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/phase',data=data_obj.data[key].ex.phase)
            exPhase.attrs['units'] = 'rad'
            exPhase.attrs['description'] = 'Phase angle of excitation force'

            exRe = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/re',data=data_obj.data[key].ex.re)
            exRe.attrs['units'] = ''
            exRe.attrs['description'] = 'Real component of excitation force'

            exIm = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/im',data=data_obj.data[key].ex.im)
            exIm.attrs['units'] = ''
            exIm.attrs['description'] = 'Imaginary component of excitation force'

            # Scattering and FK forces
            try:
                ex_sc_Mag = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/scattering/mag',data=data_obj.data[key].ex.sc.mag)
                ex_sc_Mag.attrs['units'] = ''
                ex_sc_Mag.attrs['description'] = 'Magnitude of excitation force'

                ex_sc_Phase = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/scattering/phase',data=data_obj.data[key].ex.sc.phase)
                ex_sc_Phase.attrs['units'] = 'rad'
                ex_sc_Phase.attrs['description'] = 'Phase angle of excitation force'

                ex_sc_Re = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/scattering/re',data=data_obj.data[key].ex.sc.re)
                ex_sc_Re.attrs['units'] = ''
                ex_sc_Re.attrs['description'] = 'Real component of excitation force'

                ex_sc_Im = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/scattering/im',data=data_obj.data[key].ex.sc.im)
                ex_sc_Im.attrs['units'] = ''
                ex_sc_Im.attrs['description'] = 'Imaginary component of excitation force'

                ex_fk_Mag = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/froud_krylof/mag',data=data_obj.data[key].ex.fk.mag)
                ex_fk_Mag.attrs['units'] = ''
                ex_fk_Mag.attrs['description'] = 'Magnitude of excitation force'

                ex_fk_Phase = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/froud_krylof/phase',data=data_obj.data[key].ex.fk.phase)
                ex_fk_Phase.attrs['units'] = 'rad'
                ex_fk_Phase.attrs['description'] = 'Phase angle of excitation force'

                ex_fk_Re = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/froud_krylof/re',data=data_obj.data[key].ex.fk.re)
                ex_fk_Re.attrs['units'] = ''
                ex_fk_Re.attrs['description'] = 'Real component of excitation force'

                ex_fk_Im = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/excitation/froud_krylof/im',data=data_obj.data[key].ex.fk.im)
                ex_fk_Im.attrs['units'] = ''
                ex_fk_Im.attrs['description'] = 'Imaginary component of excitation force'

            except:
                pass

            # Write added mass information
            amInf = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/inf_freq',data=data_obj.data[key].am.inf)
            amInf.attrs['units for translational degrees of freedom'] = 'kg'
            amInf.attrs['description'] = 'Infinite frequency added mass'

            am = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/all',data=data_obj.data[key].am.all)
            am.attrs['units for translational degrees of freedom'] = 'kg'
            am.attrs['units for rotational degrees of freedom'] = 'kg-m^2'
            am.attrs['description'] = 'Added mass. Frequency is the third dimension of the data structure.'

            for m in xrange(data_obj.data[key].am.all.shape[0]):

                for n in xrange(data_obj.data[key].am.all.shape[1]):

                    amComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/components/' + str(m+1) + '_' + str(n+1),data=np.array([data_obj.data[key].T, data_obj.data[key].am.all[m,n,:]]).transpose())
                    amComp.attrs['units'] = ''
                    amComp.attrs['description'] = 'Added mass components as a function of frequency'

                    radComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/components/' + str(m+1) + '_' + str(n+1),data=np.array([data_obj.data[key].T, data_obj.data[key].rd.all[m,n,:]]).transpose())
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

        rawOut = f.create_dataset('simulation_parameters/bem_raw_data',data=data_obj.data[key].bem_raw_data)
        rawOut.attrs['description'] = 'Raw output from BEM code'

        code = f.create_dataset('simulation_parameters/bem_code',data=data_obj.data[key].bem_code)
        code.attrs['description'] = 'BEM code'

        dimensional = f.create_dataset('simulation_parameters/dimensional',data=data_obj.data[key].dimensional)
        dimensional.attrs['description'] = 'True: The data is dimensional, False: The data is nondimensional'
