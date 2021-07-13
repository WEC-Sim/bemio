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


import os

import numpy as np

from bemio.data_structures import bem

try:

    from astropy.io import ascii

except:

    raise Exception('The astropy module must be installed. Try "pip install astropy"')

class NemohOutput(object):
    '''
    Class to read and interact with NEMOH simulation data

    Parameters:
        sim_dir : str, optional
            Directory where NEMOH simulation results are located
        cal_file : str, optional
            Name of NEMOH .cal file
        results_dir : float, optional
            Name of the directory that contains the NEMOH results files
        mesh_dir : str, optional
            Name of the directory that contains the NEMOH mesh files
        scale : bool, optional
            Boolean value to determine if the hydrodynamic data is scaled.
            See the bemio.data_structures.bem.scale function for more
            information

    Examples
        The following example assumes that a NEMOH simulation was run and that
        there is data ./Results and ./Mesh directories. The Nemoh.cal file is
        assumed to be located at ./Nemoh.cal

        >>> nemoh_data = NemohOtuput()
    '''
    def __init__(self, sim_dir='./', cal_file='Nemoh.cal', results_dir = 'Results', mesh_dir='Mesh', scale=False):


        # Set files
        self.scaled_at_read = scale
        self.scaled = True
        self.dir = os.path.abspath(sim_dir)

        print('\nReading NEMOH output in the ' + self.dir + ' directory')

        self.files = bem.generate_file_names(os.path.join(self.dir,cal_file))
        self.files['Nemoh']     = os.path.join(self.dir,cal_file)
        self.files['RadiationCoefficients'] = os.path.join(self.dir,results_dir,'RadiationCoefficients.tec')

        self.files['ExcitationForce'] = os.path.join(self.dir,results_dir,'ExcitationForce.tec')
        self.files['DiffractionForce'] = os.path.join(self.dir,results_dir,'DiffractionForce.tec')


        self.files['FKForce'] = os.path.join(self.dir,results_dir,'FKForce.tec')

        try:
            self.files['IRF'] = os.path.join(self.dir,results_dir,'IRF.tec')
        except:
            print('\tNo IRF forces or infinite frequency added mass forces read because the ' + self.files['IRF'] + ' was not found')

        # Initialize data ovject
        self.body = {}

        # Object to store raw data
        self.cal = bem.Raw()

        # Read cal file
        self._read_cal()

        # Read tec plot output files
        self.am, self.rd, self.w, raw_rad = _read_radiation(self.files['RadiationCoefficients'])

        self.ex_mag, self.ex_phase, temp, self.ex_mag_raw = _read_excitation(self.files['ExcitationForce'])
        self.dfr_mag, self.dfr_phase, temp, raw_diff = _read_tec(self.files['DiffractionForce'], data_type=1)
        self.fk_mag, self.fk_phase, temp, raw_fk = _read_tec(self.files['FKForce'], data_type=1)
        self.ex_phase = -1*self.ex_phase
        try:
            self.am_inf, temp1, temp2, raw_am_inf = _read_tec(self.files['IRF'], data_type=2)
        except:
            raise Exception('The file ' + self.files['IRF'] + ' was not found. Please make sure the IRF output option is set to 1 in the Nemoh.cal file')

        self.ex_im = self.ex_mag*np.sin(self.ex_phase)
        self.ex_re = self.ex_mag*np.cos(self.ex_phase)

        f_break = ['#'*100]*5
        self.raw_output = f_break + raw_rad + f_break + raw_diff + f_break + self.ex_mag_raw + f_break + raw_fk + f_break

        self._create_and_load_hydro_data_obj()

    def _create_and_load_hydro_data_obj(self):
        '''
        Function to load hydrodynamic data into HydrodynamicData object
        '''
        for i in range(self.cal.n_bods):
            self.body[i] = bem.HydrodynamicData()
            self.body[i].am.all = self.am[0+6*i:6+6*i,:]
            self.body[i].rd.all = self.rd[0+6*i:6+6*i,:]

            self.body[i].ex.mag = self.ex_mag[0+6*i:6+6*i,:,:]
            self.body[i].ex.phase = self.ex_phase[0+6*i:6+6*i,:,:]
            self.body[i].ex.im = self.ex_im[0+6*i:6+6*i,:,:]
            self.body[i].ex.re = self.ex_re[0+6*i:6+6*i,:,:]

            self.body[i].am.inf = self.am_inf[0+6*i:6+6*i, :]

            self.body[i].w = self.w
            self.body[i].T = 2.*np.pi/self.w

            self.body[i].water_depth = self.cal.water_depth
            self.body[i].g = self.cal.g
            self.body[i].rho = self.cal.rho
            self.body[i].cg = self.cal.cg[i]

            self.body[i].bem_code = 'NEMOH'
            self.body[i].bem_raw_data = self.raw_output

            self.body[i].body_num = i
            self.body[i].name = self.cal.name[i]
            self.body[i].num_bodies = self.cal.n_bods
            self.body[i].scaled = self.scaled
            self.body[i].wave_dir = self.cal.wave_dir

            self.body[i].scale(self.scaled_at_read) # Note... this is missing the KH nondimensionalization because of where it is called



    def _read_cal(self):
        '''
        Function to read Nemoh.cal file
        '''
        with open(self.files['Nemoh']) as fid:

            cal = fid.readlines()

        try:
            np.float(cal[-1].split()[0])
        except:
            cal.pop()

        self.cal.raw = cal
        self.cal.rho    = np.float(cal[1].split()[0])
        self.cal.g      = np.float(cal[2].split()[0])
        self.cal.water_depth = np.float(cal[3].split()[0])
        if self.cal.water_depth == 0:
            self.cal.water_depth = 'infinite'

        self.cal.wave_point = cal[4].split()[0:2]
        self.cal.n_bods =   int(cal[6].split()[0])

        # Read wave directions
        temp = cal[-6]
        self.cal.wave_dir_n = np.float(temp.split()[0])
        self.cal.wave_dir_start = np.float(temp.split()[1])
        self.cal.wave_dir_end = np.float(temp.split()[2])
        self.cal.wave_dir = np.linspace(self.cal.wave_dir_start,self.cal.wave_dir_end,self.cal.wave_dir_n)

        # Read frequencies
        temp = cal[-7]
        self.cal.w_n = temp.split()[0]
        self.cal.w_start = temp.split()[1]
        self.cal.w_end = temp.split()[2]

        self.cal.name = {}
        self.cal.points_panels = {}
        self.cal.n_dof = {}
        self.cal.n_forces = {}
        self.cal.n_add_lines = {}
        self.cal.dof = {}
        self.cal.forces = {}
        self.cal.add_lines = {}
        self.cal.cg = {}
        line_count = 0
        for i in range(self.cal.n_bods):
            name_with_path = cal[8+line_count].split()[0]
            self.cal.name[i] = os.path.splitext(os.path.basename(name_with_path))[0]
            self.cal.points_panels[i] = cal[9+line_count].split()[0:2]
            self.cal.n_dof[i] = int(cal[10+line_count].split()[0])
            self.cal.dof[i] = []
            self.cal.forces[i] = []
            self.cal.add_lines[i] = []
            for j in range(self.cal.n_dof[i]):
                self.cal.dof[i].append(cal[11+line_count+j])

                if int(self.cal.dof[i][-1].split()[0]) == 2:
                    self.cal.cg[i] = np.array(self.cal.dof[i][-1].split()[4:7],dtype=float)

            self.cal.n_forces[i] = int(cal[10+line_count+self.cal.n_dof[i]+1].split()[0])

            for j in range(self.cal.n_forces[i]):
                    self.cal.forces[i].append(cal[11+line_count+j+self.cal.n_dof[i]+1])

            self.cal.n_add_lines[i] = int(cal[10 + line_count + self.cal.n_dof[i] + self.cal.n_forces[i] + 2].split()[0])

            for j in range(self.cal.n_add_lines[i]):
                    self.cal.add_lines[i].append(cal[11+line_count+j+self.cal.n_dof[i]+self.cal.n_forces[i]+1])

            line_count += self.cal.n_dof[i] + self.cal.n_forces[i] + self.cal.n_add_lines[i] + 6


    def read_kh(self, file, body_num):
        '''
        Function to read NEMOH linear spring stiffness data

        Parameters:
            file : str
                Name of the file containing the linear spring stifness data
            body_num : int
                Number of the body corresponding to the Nemoh.cal input file

        Returns:
            The function does not directily return any variables, but calculates
            self.body[body_num].k (linear spring stiffness)

        Examples:
            This example assumes there is a file called `KH_1.dat`
            that contains linear spring stiffness data for body 1 in the
            `Nemoh.cal` file and that there is a NemohOutput object called
            `nemoh_data` already created.

            >>> nemoh_data.read_kh(body_num=1, file='KH_1.dat')

        .. Note:
            This function is not necessary for such a simple function, but we may
            need to make it more complicated in the future, so i'm leaving it as
            a function - mjl 25March2015
        '''
        self.body[body_num].k = np.loadtxt(file)

        if self.body[body_num].scaled is False:
            self.body[body_num].k /= (self.body[body_num].rho * self.body[body_num].g)
            print('\tSpring stiffness for body ' + self.body[body_num].name + ' scaled by read_kh method')

    def read_hydrostatics(self, file, body_num):
        '''
        Function to read NEMOH hydrostatic data

        Parameters:
            file : str
                Name of the file containing the hydrostatic data
            body_num : int
                Number of the body corresponding to the Nemoh.cal input file

        Returns:
            The function does not directily return any variables, but calculates
            self.body[body_num].disp_vol (displace volume),
            self.body[body_num].wp_area (water plane area), and
            self.body[body_num].cb (center of gravity)

        Examples:
            This example assumes there is a file called `Hydrostatics_1.dat`
            that contains hydrodynamic data for body 1 in the `Nemoh.cal` file
            and that there is a NemohOutput object called `nemoh_data` already
            created.

            >>> nemoh_data.read_hydrostatics(body_num=1,file='./Hydrostatics_1.dat')
        '''
        with open(file) as fid:

            hydrostatics = fid.readlines()

        self.body[body_num].disp_vol = np.float(hydrostatics[3].split()[-1])
        self.body[body_num].wp_area = np.float(hydrostatics[4].split()[-1])

        xf = np.float(hydrostatics[0].split()[2])
        yf = np.float(hydrostatics[1].split()[2])
        zf = np.float(hydrostatics[2].split()[2])

        self.body[body_num].cb  = np.array([xf, yf, zf])

# def _reshape_tec(data):
#     '''Internal function to reshape .tec data
#     '''
#     len = np.shape(data)[2]
#
#     out = []
#
#     for i in xrange(len):
#         out.append(data[0,:,i])
#
#     out = np.array(out)
#
#     return out

def _read_tec(file, data_type):
    '''
    Internal function to read read am and rd coefficients
    '''

    # Read added mass and damping
    with open(file) as fid:

        raw = fid.readlines()

    # Get the raw data
    proc = {}
    first = True
    for i, line in enumerate(raw):

        if 'Zone' in line:

            if first is True:
                first = False
                n_vars = i-1

            zone_length = int(line.split(',')[-2].split()[-1])
            proc[i] = ascii.read(raw[i+1:i+zone_length+1])


    # Sort the zones from the .tec file
    zones = list(proc.keys())
    zones.sort()

    # Set the frequencies and calculate number of freqs
    w = np.array(proc[zones[0]].field(0))
    n_w = np.size(w)

    if data_type == 1:
        a = np.zeros([n_vars,1,n_w])
        b = a.copy()

    if data_type == 2:
        a = np.zeros([n_vars,n_vars])
        b = []

    if data_type == 1:
        for j in range(n_vars):
            a[j,0,:] = proc[zones[-1]].field(1+j*2)
            b[j,0,:] = proc[zones[-1]].field(2+j*2)

    if data_type == 2:
        for i, zone in enumerate(zones):
            for j in range(n_vars):
                a[i,j] = proc[zone].field(1+j*2)[0]

    return (a, b, w, raw)

def _read_radiation(file, ):
    '''
    Internal function to read read am and rd coefficients
    '''

    # Read added mass and damping
    with open(file) as fid:

        raw = fid.readlines()

    # Get the raw data
    proc = {}
    first = True
    for i, line in enumerate(raw):

        if 'Zone' in line:

            if first is True:
                first = False
                n_vars = i-1

            zone_length = int(line.split(',')[-2].split()[-1])
            proc[i] = ascii.read(raw[i+1:i+zone_length+1])


    # Sort the zones from the .tec file
    zones = list(proc.keys())
    zones.sort()

    # Set the frequencies and calculate number of freqs
    w = np.array(proc[zones[0]].field(0))
    n_w = np.size(w)

    # Create and fill coefficient matrices
    a = np.zeros([n_vars,n_vars,n_w])
    b = a.copy()

    # Populate matrices
    for i, zone in enumerate(zones):
        for j in range(n_vars):
            a[i,j,:] = proc[zone].field(1+j*2)
            b[i,j,:] = proc[zone].field(2+j*2)

    return (a, b, w, raw)

def _read_excitation(file, ):
    '''
    Internal function to read read am and rd coefficients
    '''

    # Read added mass and damping
    with open(file) as fid:

        raw = fid.readlines()

    # Get the raw data
    proc = {}
    first = True
    for i, line in enumerate(raw):

        if 'Zone' in line:

            if first is True:
                first = False
                n_vars = i-1
            try:
                zone_length = int(line.split(',')[-2].split()[-1])
            except:
                pass
            try:
                zone_length = int(line.split(',')[-2].split('=')[-1])
            except:
                Exception('Tecplot format not recognized')



            proc[i] = ascii.read(raw[i+1:i+zone_length+1])


    # Sort the zones from the .tec file
    zones = list(proc.keys())
    zones.sort()

    # Set the frequencies and calculate number of freqs
    w = np.array(proc[zones[0]].field(0))
    n_w = np.size(w)

    # Create and fill coefficient matrices
    a = np.zeros([n_vars,np.size(zones),n_w])
    b = a.copy()

    # Populate matrices
    for i, zone in enumerate(zones):
        for j in range(n_vars):
            a[j,i,:] = proc[zone].field(1+j*2)
            b[j,i,:] = proc[zone].field(2+j*2)

    return (a, b, w, raw)

def read(sim_dir='./', cal_file='Nemoh.cal', results_dir = 'Results', mesh_dir='Mesh', scale=False):
    '''
    Function to read NEMOH data into a data object of type(NemohOutput)

    Parameters:
        sim_dir : str, optional
            Directory where NEMOH simulation results are located
        cal_file : str, optional
            Name of NEMOH .cal file
        results_dir : float, optional
            Name of the directory that contains the NEMOH results files
        mesh_dir : str, optional
            Name of the directory that contains the NEMOH mesh files
        scale : bool, optional
            Boolean value to determine if the hydrodynamic data is scaled.
            See the bemio.data_structures.bem.scale function for more
            information

    Returns:
        nemoh_data
            A NemohOutput object that contains the hydrodynamic data

    Examples
        The following example assumes that a NEMOH simulation was run and that
        there is data ./Results and ./Mesh directories. The Nemoh.cal file is
        assumed to be located at ./Nemoh.cal

        >>> nemoh_data = NemohOtuput()
    '''
    nemoh_data = NemohOutput(sim_dir, cal_file, results_dir, mesh_dir, scale)

    return nemoh_data

# def _reshape_tec(data):
#
#     data = data.reshape(data.shape[1],1,data.shape[0])
#
#     return out
