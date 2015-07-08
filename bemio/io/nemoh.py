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

import os

import numpy as np

from bemio.data_structures import bem

try:

    from astropy.io import ascii
 
except:

    raise Exception('The astropy module must be installed. Try "pip install astropy"')   

class NemohOutput(object):
    '''
    Class that is used to read Nemoh results.

    Inputs:
    results_dir -- the directory with the Nemoh results files (e.g. CA.dat)
    '''
    def __init__(self, sim_dir='./', cal_file='Nemoh.cal', results_dir = 'Results', mesh_dir='Mesh', out_name='nemoh_bemio.out'):

        # Set files
        self.dir = os.path.abspath(sim_dir)
        self.files = bem.generate_file_names(out_name)
        self.files['Nemoh']     = os.path.join(self.dir,cal_file)
        self.files['RadiationCoefficients'] = os.path.join(self.dir,results_dir,'RadiationCoefficients.tec')

        self.files['ExcitationForce'] = os.path.join(self.dir,results_dir,'ExcitationForce.tec')
        self.files['DiffractionForce'] = os.path.join(self.dir,results_dir,'DiffractionForce.tec')
        self.files['FKForce'] = os.path.join(self.dir,results_dir,'FKForce.tec')
        
        # Initialize data ovject
        self.data = {}

        # Object to store raw data
        self.cal = bem.Raw()

        # Read cal file
        self._read_cal()

        # Read tec plot output files
        self.am, self.rd, self.w, raw_rad = _read_tec(self.files['RadiationCoefficients'],data_type=0)

        self.ex_mag, self.ex_phase, temp, raw_ex = _read_tec(self.files['ExcitationForce'],data_type=1)
        self.dfr_mag, self.dfr_phase, temp, raw_diff = _read_tec(self.files['DiffractionForce'],data_type=1)
        self.fk_mag, self.fk_phase, temp, raw_fk = _read_tec(self.files['FKForce'],data_type=1)

        self.ex_im = self.ex_mag*np.sin(self.ex_phase)
        self.ex_re = self.ex_mag*np.cos(self.ex_phase)

        f_break = ['#'*100]*10
        self.raw_output = f_break + raw_rad + f_break + raw_diff + f_break + raw_ex + f_break + raw_fk + f_break

        self._create_and_load_hydro_data_obj()
        
    def _create_and_load_hydro_data_obj(self):
        '''
        Function to load hydrodynamic data into HydrodynamicData object
        '''
        for i in xrange(self.cal.n_bods):
            self.data[i] = bem.HydrodynamicData()
            self.data[i].am.all = self.am[0+6*i:6+6*i,:]
            self.data[i].rd.all = self.rd[0+6*i:6+6*i,:]

            self.data[i].ex.mag = self.ex_mag[:,0+6*i:6+6*i]
            self.data[i].ex.phase = self.ex_phase[:,0+6*i:6+6*i]
            self.data[i].ex.im = self.ex_im[:,0+6*i:6+6*i]
            self.data[i].ex.re = self.ex_re[:,0+6*i:6+6*i]

            self.data[i].am.inf = self.data[i].am.all[:,:,-1]

            self.data[i].w = self.w
            self.data[i].T = 2.*np.pi/self.w

            self.data[i].water_depth = self.cal.water_depth
            self.data[i].g = self.cal.g
            self.data[i].rho = self.cal.rho

            self.data[i].bem_code = 'NEMOH'
            self.data[i].bem_raw_data = self.raw_output

            self.data[i].body_num = i
            self.data[i].name = self.cal.name[i]
            self.data[i].num_bodies = self.cal.n_bods

            self.data[i].nondimensionalize_hydro_coeffs()

    def read_kh(self,body_num,file):
        '''
        Function to read HK.dat

        This function is not necessary for such a simple function, but we may
        need to make it more complicated in the future, so i'm leaving it as
        a function - mjl 25March2015
        '''
        self.data[body_num].k = np.loadtxt(file)

    def read_hydrostatics(self,body_num,file):
        '''
        Function to read hydrostatics.dat nemoh file
        '''
        with open(file) as fid:

            hydrostatics = fid.readlines()

        self.data[body_num].disp_vol = np.float(hydrostatics[3].split()[-1])
        self.data[body_num].wp_area = np.float(hydrostatics[4].split()[-1])

        xf = np.float(hydrostatics[0].split()[2])
        # xg = np.float(hydrostatics[0].split()[-1])

        yf = np.float(hydrostatics[1].split()[2])
        # yg = np.float(hydrostatics[1].split()[-1])

        zf = np.float(hydrostatics[2].split()[2])
        # zg = np.float(hydrostatics[2].split()[-1])      

        self.data[body_num].cg  = np.array([xf, yf, zf])

    def _read_cal(self):
        '''
        Function to read Nemoh.cal file
        '''
        with open(self.files['Nemoh']) as fid:

            cal = fid.readlines()

        self.cal.raw = cal
        self.cal.rho    = np.float(cal[1].split()[0])
        self.cal.g      = np.float(cal[2].split()[0])
        self.cal.n_bods =   int(cal[6].split()[0])
        self.cal.water_depth = np.float(cal[3].split()[0])

        # Read wave directions
        temp = cal[-6]
        self.cal.wave_dir_n = temp.split()[0]
        self.cal.wave_dir_start = temp.split()[1]
        self.cal.wave_dir_end = temp.split()[2]

        # Read frequencies
        temp = cal[-7]
        self.cal.w_n = temp.split()[0]
        self.cal.w_start = temp.split()[1]
        self.cal.w_end = temp.split()[2]

        add_lines = 0
        self.cal.name = {}
        for i in xrange(self.cal.n_bods):
            self.cal.name[i] = cal[8+i*18+add_lines].split()[0]
            add_lines += int(cal[24+18*i+add_lines].split()[0])


            


def _reshape_tec(data):

    len = np.shape(data)[2]

    out = []

    for i in xrange(len):
        out.append(data[0,:,i])

    out = np.array(out)

    return out

def _read_tec(file,data_type=0):
    '''
    Function to read read am and rd coefficients

    Internal function called during at __init__

    Need to describe data_type
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
    zones = proc.keys()
    zones.sort()

    # Set the frequencies and calculate number of freqs
    w = np.array(proc[zones[0]].field(0))
    n_w = np.size(w)

    # Create and fill coefficient matrices
    a = np.zeros([n_vars,n_vars,n_w])
    b = a.copy()

    # Populate matrices
    for i, zone in enumerate(zones):

        for j in xrange(n_vars):
            a[i,j,:] = proc[zone].field(1+j*2)     
            b[i,j,:] = proc[zone].field(2+j*2)

    if data_type == 1:
        a = _reshape_tec(a)
        b = _reshape_tec(b)

    return (a, b, w, raw)

