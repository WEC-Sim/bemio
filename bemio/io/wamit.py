# Copyright 2014 the National Renewable Energy Laboratory and Sandia Corporation
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
'''bemio WAMIT module

This moduel provides functionality to read and interact with WAMIT simulation
output data
'''
import os

import numpy as np

from bemio.data_structures import bem

class WamitOutput(object):
    '''
    Class to read and interact with WAMIT simulation data

    Parameters:
        out_file : str
            Name of the wamit .out output file. In order to read scattering
            and Froud-Krylof forces the .3sc (scatterin) and .3fk
            (Froud-Krylof) coefficinet files must have the same base name
            as the .out file.
        density : float, optional
            Water density used to scale the hydrodynamic coefficient data
        gravity : float, optional
            Acceleration due to gravity used to scale the hydrodynamic
            coefficient dataaaa
        scale : bool, optional
            Boolean value to determine if the hydrodynamic data is scaled.
            See the bemio.data_structures.bem.scale function for more
            information
        ex_calc : str
            Excitation coefficients to be used, 'diffraction' or 'haskind'.

    Examples
        The user can create a WamitOutput data object directily as show below,
        or the bemio.io.wamit.read function can be used. The following example
        assumes there is a WAMIT output file named `wamit.out`.

        >>> wamit_data = WamitOtuput(out_file=wamit.out)
    '''
    def __init__(self, out_file, density=1000., gravity=9.81, scale=False, ex_calc='diffraction'):

        self.files = bem.generate_file_names(out_file)
        self.files['3sc'] = self.files['base_name'] + '.3sc'
        self.files['3fk'] = self.files['base_name'] + '.3fk'

        self.rho = density
        self.g = gravity

        self.body = {}
        self.scaled_at_read = scale
        self.scaled = False
        self.ex_calc = ex_calc
        self._read()

    def _read(self):
        '''Internal function to read WAMIT output file into the class. that is called during __init__
        '''

        print('\nReading the WAMIT results in the ' + self.files['out'] + ' file')

        with open(self.files['out'],'rU') as fid:

            raw = fid.readlines()

        code = 'WAMIT'
        num_bodies = 0 # Total number of bodies
        bod_count = 0 # Counter for bodies
        T = []
        cg = {}
        cb = {}
        name = {}
        disp_vol = {}
        k = {}
        wave_dir = []
        empty_line = '\n'


        for i, line in enumerate(raw):


            if "POTEN run date and starting time:" in line:

                skip = 2
                data = raw[i+skip]
                count = 0

                while data != empty_line:

                    if float(data.split()[0]) == 0. or float(data.split()[0]) == -1.:
                        count += 1
                        data = raw[i+count+skip]

                    else:
                        count += 1
                        T.append(float(data.split()[0]))
                        data = raw[i+count+skip]

            if "Wave Heading (deg)" in line:
                wave_dir.append(float(line.split()[-1]))

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

                for j in range(20): # look for bodies within the next 20 lines

                    if "N=" in raw[i+j]:

                        num_bodies += 1
                        name[num_bodies-1] = raw[i+j].split()[-1]


            # Read the body positions
            # if "Total panels:" in line or "NPATCH:" in line:
            if "XBODY" in line:
                for j in range(12): # look for position within the next 15 lines - will only work for wamit files of about 5 bodies

                    if 'XBODY =' in raw[i+j]:
                        '''
                        Note that this is the XBOD YBOD ZBOD defined in the wamit .out file, not the cg as defined in the wamit file
                        '''

                        temp = raw[i+j].split()
                        cg[bod_count] = np.array([temp[2],temp[5],temp[8]]).astype(float)

                    if 'Volumes (VOLX,VOLY,VOLZ):' in raw[i+j]:

                        temp = raw[i+j].split()
                        disp_vol[bod_count] = float(temp[-1])

                    if 'Center of Buoyancy (Xb,Yb,Zb):' in raw[i+j]:

                        temp = raw[i+j].split()
                        cb[bod_count] = np.array([temp[-3],temp[-2],temp[-1]]).astype(float)

                    if 'C(3,3),C(3,4),C(3,5):' in raw[i+j]:

                        temp = np.zeros([6,6])
                        temp2 = raw[i+j].split()
                        temp[2,2] = np.float(temp2[1])
                        temp[2,3] = np.float(temp2[2])
                        temp[2,4] = np.float(temp2[3])
                        temp[3,2] = temp[2,3]
                        temp[4,2] = temp[2,4]

                        temp2 = raw[i+j+1].split()
                        temp[3,3] = np.float(temp2[1])
                        temp[3,4] = np.float(temp2[2])
                        temp[3,5] = np.float(temp2[3])
                        temp[4,3] = temp[3,4]
                        temp[5,3] = temp[3,5]

                        temp2 = raw[i+j+2].split()
                        temp[4,4] = np.float(temp2[1])
                        temp[4,5] = np.float(temp2[2])
                        temp[5,4] = temp[4,5]

                        k[bod_count] = temp


                bod_count += 1

        # Put things into numpy arrays
        T = np.array(T).astype(float)
        wave_dir = np.array(wave_dir).astype(float)

        # Only select the wave headings once
        temp = 999999
        temp_wave_dir = []
        count = 0

        while temp != wave_dir[0]:

            count += 1
            temp_wave_dir.append(wave_dir[count-1])
            temp = wave_dir[count]


        wave_dir = np.array(temp_wave_dir).astype(float)

        # Read added mass and rad damping
        count_freq = 0
        am_all = np.zeros([6*num_bodies,6*num_bodies,T.size])
        rd_all = am_all.copy()
        am_inf = np.zeros([6*num_bodies,6*num_bodies])
        am_zero = am_inf.copy()

        for i, line in enumerate(raw):

            # Read inf freq added mass
            if "Wave period = zero" in line:

                count = 7
                temp_line = raw[count+i]

                while temp_line != empty_line:

                    am_inf[int(temp_line.split()[0])-1,int(temp_line.split()[1])-1] = temp_line.split()[2]
                    count += 1
                    temp_line = raw[count+i]


            # Read zero freq added mass
            if "Wave period = infinite" in line:

                count = 7
                temp_line = raw[count+i]

                while temp_line != empty_line:

                    am_zero[int(temp_line.split()[0])-1,int(temp_line.split()[1])-1] = temp_line.split()[2]
                    count += 1
                    temp_line = raw[count+i]


            # Read freq dependent added mass and rad damping
            if "Wave period (sec) =" in line:
                temp = raw[i].split()
                T[count_freq]=temp[4]

                count = 7
                temp_line = raw[count+i]

                while temp_line != empty_line:

                    am_all[int(temp_line.split()[0])-1,int(temp_line.split()[1])-1,count_freq] = temp_line.split()[2]
                    rd_all[int(temp_line.split()[0])-1,int(temp_line.split()[1])-1,count_freq] = temp_line.split()[3]
                    count += 1
                    temp_line = raw[count+i]

                count_freq += 1


        # Terribly complicated code to read excitation forces and phases, RAOs, etc
        ex_all = np.zeros([6*num_bodies,wave_dir.size,T.size])
        phase_all = ex_all.copy()
        rao_all = ex_all.copy()
        rao_phase_all = ex_all.copy()
        ssy_all = ex_all.copy()
        ssy_phase_all = ex_all.copy()
        haskind_all = ex_all.copy()
        haskind_phase_all = ex_all.copy()
        count_diff2 = 0
        count_rao2 = 0
        count_ssy2 = 0
        count_haskind2 = 0
        for i, line in enumerate(raw):

            count_diff = 0
            count_rao = 0
            count_ssy = 0
            count_haskind = 0

            if "DIFFRACTION EXCITING FORCES AND MOMENTS" in line:

                count_diff += 1
                count_diff2 += 1
                count_wave_dir = 0
                count = 0

                while count_wave_dir < wave_dir.size:

                    count += 1

                    if "Wave Heading (deg) :" in raw[i+count_diff + count]:

                        count_wave_dir += 1
                        temp_line = raw[i+count_diff+count+4]
                        count2 = 0

                        while temp_line != empty_line:
                            count2 += 1
                            ex_all[int(temp_line.split()[0])-1,count_wave_dir-1,count_diff2-1] = float(temp_line.split()[1])
                            phase_all[int(temp_line.split()[0])-1,count_wave_dir-1,count_diff2-1] = float(temp_line.split()[2])
                            temp_line = raw[i+count_diff+count+4+count2]

            if "RESPONSE AMPLITUDE OPERATORS" in line:

                count_rao += 1
                count_rao2 += 1
                count_wave_dir = 0
                count = 0

                while count_wave_dir < wave_dir.size:

                    count += 1

                    if "Wave Heading (deg) :" in raw[i+count_rao + count]:

                        count_wave_dir += 1
                        temp_line = raw[i+count_rao+count+4]
                        count2 = 0

                        while temp_line != empty_line:
                            count2 += 1
                            rao_all[int(temp_line.split()[0])-1,count_wave_dir-1,count_rao2-1] = float(temp_line.split()[1])
                            rao_phase_all[int(temp_line.split()[0])-1,count_wave_dir-1,count_rao2-1] = float(temp_line.split()[2])
                            temp_line = raw[i+count_rao+count+4+count2]


            if "HASKIND EXCITING FORCES AND MOMENTS" in line:

                count_haskind += 1
                count_haskind2 += 1
                count_wave_dir = 0
                count = 0

                while count_wave_dir < wave_dir.size:

                    count += 1

                    if "Wave Heading (deg) :" in raw[i+count_haskind + count]:

                        count_wave_dir += 1
                        temp_line = raw[i+count_haskind+count+4]
                        count2 = 0

                        while temp_line != empty_line:
                            count2 += 1
                            haskind_all[int(temp_line.split()[0])-1,count_wave_dir-1,count_haskind2-1] = float(temp_line.split()[1])
                            haskind_phase_all[int(temp_line.split()[0])-1,count_wave_dir-1,count_haskind2-1] = float(temp_line.split()[2])
                            temp_line = raw[i+count_haskind+count+4+count2]

            if "SURGE, SWAY & YAW DRIFT FORCES (Momentum Conservation)" in line:

                count_ssy += 1
                count_ssy2 += 1
                count_wave_dir = 0
                count = 0

                while count_wave_dir < wave_dir.size:

                    count += 1

                    if "Wave Heading (deg) :" in raw[i+count_ssy + count]:

                        count_wave_dir += 1
                        temp_line = raw[i+count_ssy+count+4]
                        count2 = 0

                        while temp_line != empty_line:
                            count2 += 1
                            ssy_all[int(temp_line.split()[0])-1,count_wave_dir-1,count_ssy2-1] = float(temp_line.split()[1])
                            ssy_phase_all[int(temp_line.split()[0])-1,count_wave_dir-1,count_ssy2-1] = float(temp_line.split()[2])
                            temp_line = raw[i+count_ssy+count+4+count2]

        if os.path.exists(self.files['3sc']):
            sc_re = np.zeros([6*num_bodies,wave_dir.size,T.size])
            sc_im = sc_re.copy()
            sc_phase = sc_re.copy()
            sc_mag = sc_re.copy()

            scattering = np.loadtxt(self.files['3sc'],skiprows=1)
            line_count = 0

            for freq_n in range(T.size):
                for beta_n in range(wave_dir.size):
                    wave_dir_hold =  scattering[line_count][1]
                    while line_count < scattering.shape[0] and scattering[line_count][1] == wave_dir_hold:
                        comp = int(scattering[line_count][2])-1
                        sc_mag[comp,beta_n,freq_n] = scattering[line_count][3]
                        sc_phase[comp,beta_n,freq_n] = scattering[line_count][4]
                        sc_re[comp,beta_n,freq_n] = scattering[line_count][5]
                        sc_im[comp,beta_n,freq_n] = scattering[line_count][6]
                        wave_dir_hold =  scattering[line_count][1]
                        line_count += 1

        else:
            print('\tThe file ' + self.files['3sc'] + ' does not exist... not reading scattering coefficients.')

        if os.path.exists(self.files['3fk']):
            fk_re = np.zeros([6*num_bodies,wave_dir.size,T.size])
            fk_im = fk_re.copy()
            fk_phase = fk_re.copy()
            fk_mag = fk_re.copy()

            fk = np.loadtxt(self.files['3fk'],skiprows=1)
            line_count = 0

            for freq_n in range(T.size):
                for beta_n in range(wave_dir.size):
                    wave_dir_hold =  fk[line_count][1]
                    while line_count < fk.shape[0] and fk[line_count][1] == wave_dir_hold:
                        comp = int(fk[line_count][2])-1
                        fk_mag[comp,beta_n,freq_n] = fk[line_count][3]
                        fk_phase[comp,beta_n,freq_n] = fk[line_count][4]
                        fk_re[comp,beta_n,freq_n] = fk[line_count][5]
                        fk_im[comp,beta_n,freq_n] = fk[line_count][6]
                        wave_dir_hold =  scattering[line_count][1]
                        line_count += 1

        else:
            print('\tThe file ' + self.files['3fk'] + ' does not exist... not reading froud krylof coefficients.')

        # Load data into the hydrodata structure
        for i in range(num_bodies):
            self.body[i] = bem.HydrodynamicData()
            self.body[i].scaled = self.scaled
            self.body[i].g = self.g
            self.body[i].rho = self.rho
            self.body[i].body_num = i
            self.body[i].name = name[i][0:-4]
            self.body[i].water_depth = water_depth
            self.body[i].num_bodies = num_bodies
            self.body[i].cg = cg[i]
            self.body[i].cb = cb[i]
            self.body[i].k = k[i]
            self.body[i].disp_vol = disp_vol[i]
            self.body[i].wave_dir = wave_dir
            self.body[i].T = T
            self.body[i].w = 2.0*np.pi/self.body[i].T

            if 'am_inf' in locals():

                self.body[i].am.inf = am_inf[6*i:6+6*i,:]

            else:

                self.body[i].am.inf = np.nan*np.zeros([6*num_bodies,6*num_bodies,self.body[i].T.size])
                print('Warning: body ' + str(i) + ' - The WAMTI .out file specified does not contain infinite frequency added mass coefficients')


            if 'am_zero' in locals():

                self.body[i].am.zero = am_zero[6*i:6+6*i,:]

            else:

                self.body[i].am.zero = np.nan*np.zeros([6*num_bodies,6*num_bodies,self.body[i].T.size])
                print('Warning: body ' + str(i) + ' - The WAMTI .out file specified does not contain zero frequency added mass coefficients')


            if 'am_all' in locals():

                self.body[i].am.all = am_all[6*i:6+6*i,:,:]
            else:

                self.body[i].am.all = np.nan*np.zeros([6*num_bodies,6*num_bodies,self.body[i].T.size])
                print('Warning: body ' + str(i) + ' - The WAMTI .out file specified does not contain any frequency dependent added mass coefficients')


            if 'rd_all' in locals():

                self.body[i].rd.all = rd_all[6*i:6+6*i,:,:]

            else:

                self.body[i].rd.all = np.nan*np.zeros([6*num_bodies,6*num_bodies,self.body[i].T.size])
                print('Warning: body ' + str(i) + ' - The WAMTI .out file specified does not contain any frequency dependent radiation damping coefficients')

            if 'ex_all' in locals() and self.ex_calc == 'diffraction' :
                
                self.body[i].ex.mag = ex_all[6*i:6+6*i,:,:]
                self.body[i].ex.phase = np.deg2rad(phase_all[6*i:6+6*i,:,:])
                self.body[i].ex.re = self.body[i].ex.mag*np.cos(self.body[i].ex.phase)
                self.body[i].ex.im = self.body[i].ex.mag*np.sin(self.body[i].ex.phase)
            
            elif 'haskind_all' in locals() and self.ex_calc == 'haskind' :
                
                self.body[i].ex.mag = haskind_all[6*i:6+6*i,:,:]
                self.body[i].ex.phase = np.deg2rad(haskind_phase_all[6*i:6+6*i,:,:])
                self.body[i].ex.re = self.body[i].ex.mag*np.cos(self.body[i].ex.phase)
                self.body[i].ex.im = self.body[i].ex.mag*np.sin(self.body[i].ex.phase)

            else:

                print('Warning: body ' + str(i) + ' - The WAMTI .out file specified does not contain any excitation coefficients')

            if 'sc_mag' in locals():
                self.body[i].ex.sc.mag = sc_mag[6*i:6+6*i,:,:]
                self.body[i].ex.sc.phase = np.deg2rad(sc_phase[6*i:6+6*i,:,:])
                self.body[i].ex.sc.re = sc_re[6*i:6+6*i,:,:]
                self.body[i].ex.sc.im = sc_im[6*i:6+6*i,:,:]

            else:
                pass
                # print 'Warning: body ' + str(i) + ' - The WAMTI .3sc file specified does not contain any scattering coefficients'

            if 'fk_mag' in locals():
                self.body[i].ex.fk.mag = fk_mag[6*i:6+6*i,:,:]
                self.body[i].ex.fk.phase = np.deg2rad(fk_phase[6*i:6+6*i,:,:])
                self.body[i].ex.fk.re = fk_re[6*i:6+6*i,:,:]
                self.body[i].ex.fk.im = fk_im[6*i:6+6*i,:,:]

            else:
                pass
                # print 'Warning: body ' + str(i) + ' - The WAMTI .3fk file specified does not contain any froude krylof coefficients'

            if 'rao_all' in locals():

                self.body[i].rao.mag = rao_all[6*i:6+6*i,:,:]
                self.body[i].rao.phase = np.deg2rad(phase_all[6*i:6+6*i,:,:])
                self.body[i].rao.re = self.body[i].rao.mag*np.cos(self.body[i].rao.phase)
                self.body[i].rao.im = self.body[i].rao.mag*np.sin(self.body[i].rao.phase)

            else:

                print('Warning: body ' + str(i) + ' - The WAMTI .out file specified does not contain any rao data')

            if 'ssy_all' in locals():

                self.body[i].ssy.mag = ssy_all[6*i:6+6*i,:,:]
                self.body[i].ssy.phase = np.deg2rad(phase_all[6*i:6+6*i,:,:])
                self.body[i].ssy.re = self.body[i].ssy.mag*np.cos(self.body[i].ssy.phase)
                self.body[i].ssy.im = self.body[i].ssy.mag*np.sin(self.body[i].ssy.phase)

            else:

                print('Warning: body ' + str(i) + ' - The WAMTI .out file specified does not contain any rao data')

            self.body[i].bem_raw_data = raw
            self.body[i].bem_code = code

            self.body[i].scale(scale=self.scaled_at_read)

def read(out_file, density=1000., gravity=9.81, scale=False, ex_calc='diffraction'):
    '''
    Function to read WAMIT data into a data object of type(WamitOutput)

    Parameters:
        out_file : str
            Name of the wamit .out output file. In order to read scattering
            and Froud-Krylof forces the .3sc (scatterin) and .3fk
            (Froud-Krylof) coefficinet files must have the same base name
            as the .out file.
        density : float, optional
            Water density used to scale the hydrodynamic coefficient data
        gravity : float, optional
            Acceleration due to gravity used to scale the hydrodynamic
            coefficient data
        scale : bool, optional
            Boolean value to determine if the hydrodynamic data is scaled.
            See the bemio.data_structures.bem.scale function for more
            information
        ex_calc : str
            Excitation coefficients to be used, 'diffraction' or 'haskind'.

    Returns:
        wamit_data
            A WamitData object that contains the data from the WAMIT .out
            file specified

    Examples:
        The following example assumes there is a WAMIT output file named
        `wamit.out`

        >>> wamit_data = read(out_file=wamit.out)
    '''
    wamit_data = WamitOutput(out_file, density, gravity, scale, ex_calc)

    return wamit_data
