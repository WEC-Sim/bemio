# Copyright 2014 the National Renewable Energy Laboratory
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

"""
This module serves the purpose of reading and writing  surface meshes from the
following mesh formats:
    * WAMTI
    * Nemoh
    * VTK Polydata (*.vtp)
    * STL files

Limitations: Currentoy the writing functions assume that all meshes are quad
meshes. A few simple improvements (i.e. if tests) to the writing functions can fix this.

Author: Michael Lawson
"""
import numpy as np

import os

from sys import platform as _platform

try:

    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

except:

    raise Exception('The VTK Python module is required to use this module.')


class PanelMesh(object):
    ''' Class to handle surface meshes. The class currently handles quad or tri
    meshes
    '''
    def __init__(self,file_name):

        self.file_name = file_name
        self.faces = []
        self.points = []
        self.num_faces = None
        self.num_points = None
        self.orig_type = None
        self._center_of_gravity = None
        self._center_of_buoyancy = None
        self._volume = None
        self._volume_x = None
        self._volume_y = None
        self._volume_z = None
        self._surface_area = None
        self._normals = None
        self._cell_surface_area = None
        self._centroid = None
        self._hydrostatic_stiffness = None

        self.water_density = 1000.

        if os.path.isfile(file_name) is False:

            raise Exception('The file ' + file_name + ' does not exist')

    def __repr__(self):

        out_string = 'File name: ' + str(self.file_name) + \
        '\nNumber of points: ' + str(self.num_points) + \
        '\nNumber of faces: ' + str(self.num_faces) + \
        '\nOriginal mesh type: ' + str(self.orig_type) + \
        '\nObject type: bemio.mesh_utilities.mesh.PanelMesh'

        return out_string

    # @property
    # def center_of_gravity(self, ):
    #     if self._center_of_mass is None:
    #
    #         # com = vtk.vtkCenterOfMass()
    #         # com.SetInputData(self.vtp_mesh)
    #         # com.Update()
    #         # self._center_of_mass = com.GetCenter()
    #
    #     return self._center_of_mass

    @property
    def hydrostatic_stiffness(self, ):
        if self._hydrostatic_stiffness is None:

            self._hydrostatic_stiffness = np.zeros([6,6])

            for face_n,face in enumerate(self.faces):
                if np.sum(self.points[face[0]][2]+self.points[face[1]][2]+self.points[face[2]][2]+self.points[face[3]][2]) != 0.:
                    self._hydrostatic_stiffness[2,2] += -self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[2,3] += -self.centroid[face_n][1] * self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[2,4] += self.centroid[face_n][0] * self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[3,3] += -self.centroid[face_n][1]**2 * self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[3,4] += self.centroid[face_n][0] * self.centroid[face_n][1] * self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[4,4] += -self.centroid[face_n][0]**2 * self.normals[face_n][2] * self.cell_surface_area[face_n]

            self._hydrostatic_stiffness[3,3] += self.volume * self.center_of_buoyancy[2] - self.volume * self.center_of_gravity[2]
            self._hydrostatic_stiffness[3,5] += -self.volume * self.center_of_buoyancy[0] + self.volume * self.center_of_gravity[0]
            self._hydrostatic_stiffness[4,4] += self.volume * self.center_of_buoyancy[2] - self.volume * self.center_of_gravity[2]
            self._hydrostatic_stiffness[4,5] += -self.volume * self.center_of_buoyancy[1] + self.volume * self.center_of_gravity[1]

        return self._hydrostatic_stiffness

    @property
    def center_of_buoyancy(self, ):
        if self._center_of_buoyancy is None:
            x_b = 0.
            y_b = 0.
            z_b = 0.
            self._center_of_buoyancy = 0.

            for face_n in xrange(self.num_faces):
                x_b += self.normals[face_n][0]*self.centroid[face_n][0]**2*self.cell_surface_area[face_n]
                y_b += self.normals[face_n][1]*self.centroid[face_n][1]**2*self.cell_surface_area[face_n]
                z_b += self.normals[face_n][2]*self.centroid[face_n][2]**2*self.cell_surface_area[face_n]
            self._center_of_buoyancy = 1./(2.*self.volume)*np.array([x_b, y_b, z_b])

        return self._center_of_buoyancy

    @property
    def normals(self, ):
        if self._normals is None:
            self._normals = {}
            for face_n in xrange(self.num_faces):
                a = self.points[self.faces[face_n][1]] - self.points[self.faces[face_n][0]]
                b = self.points[self.faces[face_n][2]] - self.points[self.faces[face_n][1]]
                self._normals[face_n] = np.cross(a,b)

                if self._normals[face_n][0] == 0. and  self._normals[face_n][1] == 0. and self._normals[face_n][2] == 0.:
                    a = self.points[self.faces[face_n][2]] - self.points[self.faces[face_n][1]]
                    b = self.points[self.faces[face_n][3]] - self.points[self.faces[face_n][2]]
                    self._normals[face_n] = np.cross(a,b)

                if self._normals[face_n][0] == 0. and  self._normals[face_n][1] == 0. and self._normals[face_n][2] == 0.:
                    a = self.points[self.faces[face_n][2]] - self.points[self.faces[face_n][0]]
                    b = self.points[self.faces[face_n][3]] - self.points[self.faces[face_n][2]]
                    self._normals[face_n] = np.cross(a,b)

                self._normals[face_n] /= np.linalg.norm(self._normals[face_n])

        return self._normals

    @property
    def cell_surface_area(self):
        if self._cell_surface_area is None:
            self._cell_surface_area = {}
            for face_n in xrange(self.num_faces):
                a = self.points[self.faces[face_n][1]] - self.points[self.faces[face_n][0]]
                b = self.points[self.faces[face_n][2]] - self.points[self.faces[face_n][1]]
                c = self.points[self.faces[face_n][3]] - self.points[self.faces[face_n][2]]
                d = self.points[self.faces[face_n][0]] - self.points[self.faces[face_n][3]]
                self._cell_surface_area[face_n] = 1./2. * ( np.linalg.norm(np.cross(a,b)) + np.linalg.norm(np.cross(c,d)) )

        return self._cell_surface_area

    @property
    def volume(self):
        if self._volume is None:
            tri_converter = vtk.vtkTriangleFilter()
            tri_converter.SetInputDataObject(self.vtp_mesh)
            tri_converter.Update()
            tri_mesh = tri_converter.GetOutput()
            mass_props = vtk.vtkMassProperties()
            mass_props.SetInputDataObject(tri_mesh)
            self._volume = mass_props.GetVolume()

        return self._volume

    @property
    def centroid(self):

        if self._centroid is None:
            self._centroid = {}

            for face_n,face in enumerate(self.faces):
                points = [self.points[face[0]],
                          self.points[face[1]],
                          self.points[face[2]],
                          self.points[face[3]]]

                points = map(np.asarray, set(map(tuple, points))) # This removes duplicate points
                self._centroid[face_n] = np.mean(points,axis=0)

        return self._centroid

    @property
    def volume_x(self):
        if self._volume_x is None:
            self._calc_component_vol()
        return self._volume_x

    @property
    def volume_y(self):
        if self._volume_y is None:
            self._calc_component_vol()
        return self._volume_y

    @property
    def volume_z(self):
        if self._volume_z is None:
            self._calc_component_vol()
        return self._volume_z

    def _calc_component_vol(self, ):
        self._volume_x = 0.
        self._volume_y = 0.
        self._volume_z = 0.
        volume = 0.

        for face_n in xrange(self.num_faces):
            volume += self.normals[face_n]*self.centroid[face_n]*self.cell_surface_area[face_n]
            # self._volume_x += self.normals[face_n][0]*self.centroid[face_n][0]*self.cell_surface_area[face_n]
            # self._volume_y += self.normals[face_n][1]*self.centroid[face_n][1]*self.cell_surface_area[face_n]
            # self._volume_z += self.normals[face_n][2]*self.centroid[face_n][2]*self.cell_surface_area[face_n]

        self._volume_x = volume[0]
        self._volume_y = volume[1]
        self._volume_z = volume[2]

    @property
    def surface_area(self):
        tri_converter = vtk.vtkTriangleFilter()
        tri_converter.SetInputDataObject(self.vtp_mesh)
        tri_converter.Update()
        tri_mesh = tri_converter.GetOutput()
        mass_props = vtk.vtkMassProperties()
        mass_props.SetInputDataObject(tri_mesh)
        self._surface_area = mass_props.GetSurfaceArea()
        return self._surface_area

    def _create_vtp_mesh(self):
        self.vtp_mesh    = vtk.vtkPolyData()
        points  = vtk.vtkPoints()
        polys   = vtk.vtkCellArray()

        for i in range(self.num_points):

            points.InsertPoint(i, self.points[i])


        for i in range(self.num_faces):

            polys.InsertNextCell(_mk_vtk_id_list(self.faces[i]))

        self.vtp_mesh.SetPoints(points)
        self.vtp_mesh.SetPolys(polys)

    def _collapse(self,plane=2,value=0.0,direction=1):
        '''Collapse points
        '''
        for i in xrange(self.num_faces):

            for j in xrange(self.faces[i].size):

                p = int(self.faces[i][j])

                if self.points[p][plane] > value*direction:

                    self.points[p][plane] = value

    def write(self,mesh_format='VTK'):

        out_file_base = os.path.splitext(self.file_name)[0] + '_bemio_output'

        if mesh_format == 'VTK' or mesh_format == 'VTP':
            self.out_file = out_file_base + '.vtp'
            self._write_vtp()

        if mesh_format == 'WAMIT' or mesh_format == 'GDF':
            self.out_file = out_file_base + '.gdf'
            self._write_gdf()

        if mesh_format == 'NEMOH':
            self.out_file = out_file_base + '.dat'
            self._write_nemoh()


    def _write_vtp(self):

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(self.out_file)
        writer.SetInputData(self.vtp_mesh)
        writer.SetDataModeToAscii()
        writer.Write()

        print 'Wrote VTK PolyData mesh to ' + str(self.out_file)

    def _write_nemoh(self):

        with open(self.out_file,'w') as fid:
            fid.write('2 0') # This should not be hard coded
            fid.write('\n')
            for i in xrange(self.num_points):
                fid.write(str(i+1) + ' ' +str(self.points[i]).replace('[','').replace(']',''))
                fid.write('\n')
            fid.write('0 0 0 0')
            fid.write('\n')
            for i in xrange(self.num_faces):
                fid.write(str(self.faces[i]+1).replace('[','').replace(']','').replace('.',''))
                fid.write('\n')
            fid.write('0 0 0 0')

        print 'Wrote NEMOH mesh to ' + str(self.out_file)


    def _write_gdf(self,out_file=None):

        with open(self.out_file,'w') as fid:
            fid.write('Mesh file written by meshio.py')
            fid.write('\n')
            fid.write('1 9.80665       ULEN GRAV')
            fid.write('\n')
            fid.write('0  0    ISX  ISY')
            fid.write('\n')
            fid.write(str(self.num_faces))
            fid.write('\n')
            for i,face in enumerate(self.faces):
                if np.size(face) is 4: # if the mesh element is a quad
                    for j,pointKey in enumerate(face):
                        fid.write(str(self.points[pointKey]).replace(',','').replace('[','').replace(']','') + '\n')
                if np.size(face) is 3: # if the mesh element is a tri
                    faceMod = np.append(face,face[-1])
                    for j,pointKey in enumerate(faceMod):
                        fid.write(str(self.points[pointKey]).replace(',','').replace('[','').replace(']','') + '\n')

        print 'Wrote WAMIT mesh to ' + str(self.out_file)


    def cut(self,plane=2,value=0.0,direction=1):
        self.collapse(plane,value,direction)

        tempFaces = []
        count = 0

        for i in xrange(self.num_faces):

           delete_face = 0

           for j in xrange(4):

               p = self.faces[i][j]
               z = float(self.cords[int(p)][2])

               if z == 0.:
                   delete_face += 1

           if delete_face != 4:
               tempFaces.append(self.faces[i])
               count  += 1

        print 'removed ' + str(count) + ' surface faces'
        self.faces = tempFaces
        self.num_faces = self.faces.shape[0]

    def view(self,color=[0.5,1,0.5],opacity=1.0):

        # Create a mapper and load VTP data into the mapper
        mapper=vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.vtp_mesh)

        # Create an actor that contains the data in the mapper
        actor=vtk.vtkActor()
        actor.GetProperty().SetColor([0.5,1,0.5])
        actor.GetProperty().SetOpacity(1.0)
        actor.SetMapper(mapper)
        actor.GetProperty().EdgeVisibilityOn()

        # Add axes
        axes = vtk.vtkAxesActor()

        # Render the data
        ren = vtk.vtkRenderer()
        ren.AddActor(actor)
        ren.AddActor(axes)

        # Create a render window
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        renWin.SetSize(600, 600)

        # Start the visiuilization
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        ren.SetBackground(0,0,0)
        renWin.Render()
        iren.Start()

        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

    def paraview(self):
        '''
        Visualize the geometry in paraview - this is kind of a hack function

        To use this function make a symbolic link to the paraview.app folder
        on your system to ~/bin. Or alternatively change this function
        '''

        file_name = self.file_name[:-3] + 'vis-temp.vtp'
        if _platform == 'darwin':


            self.write_vtp(out_file=file_name)
            try:

                os.system('open ~/bin/paraview ' + file_name + ' &')
            except:

                raise Exception('~/bin/paraview not found')
        else:

            print 'paraview() function only supported for osx'


def _read_gdf(file_name):
    '''
    Function to read WAMIT GDF meshes

    Inputs:
        file_name: name of the mesh file to be read
    Outputs:
        mesh_data: panelMesh object containing the mesh data
    '''

    with open(file_name,'r') as fid:

        lines = fid.readlines()

    mesh_data = PanelMesh(file_name)
    mesh_data.orig_type = 'WAMIT (.gdf)'

    mesh_data.gdfLines = lines
    mesh_data.uLen = int(lines[1].split()[0])
    mesh_data.gravity = float(lines[1].split()[1])
    mesh_data.isx = float(lines[2].split()[0])
    mesh_data.isy = float(lines[2].split()[1])
    mesh_data.num_faces = int(lines[3].split()[0])
    mesh_data.num_points = mesh_data.num_faces * 4
    mesh_data.points = np.array([temp.split() for temp in lines[4:]]).astype(np.float)

    mesh_data.pointsString = [str(temp).replace(","      ,'').replace('\r','') for temp in lines[4:]] # Output string for Nemoh mesh fil

    for panelNum,i in enumerate(np.arange(4,4+mesh_data.num_points,4)):

        mesh_data.faces.append(np.array([i-4,i-3,i-2,i-1]))

    print 'Successfully read WAMIT mesh ' + str(file_name)

    return mesh_data


def _read_stl(file_name):
    '''
    Function to read STL meshes

    Inputs:
        file_name: name of the mesh file to be read
    Outputs:
        mesh_data: panelMesh object containing the mesh data
    '''

    reader = vtk.vtkSTLReader()
    reader.SetFileName(file_name)
    reader.Update()

    mesh_data = PanelMesh(file_name)
    mesh_data.orig_type = 'Stereolithography (.stl)'
    mesh_data.num_faces = int(reader.GetOutput().GetNumberOfCells())
    mesh_data.num_points = mesh_data.num_faces * 3
    for i in range(mesh_data.num_faces):
        n = i*3
        mesh_data.faces.append(np.array([n,n+1,n+2,n+2]))
        mesh_data.points.append(np.array(vtk_to_numpy(reader.GetOutput().GetCell(i).GetPoints().GetData())))
    mesh_data.points = np.array(mesh_data.points).reshape([mesh_data.num_faces*3,3])


    print 'Successfully read STL mesh ' + str(file_name)


    return mesh_data


def _read_vtp(file_name):
    '''
    Function to read VTK Polydata meshes

    Inputs:
        file_name: name of the mesh file to be read
    Outputs:
        mesh_data: panelMesh object containing the mesh data
    '''

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_name)
    reader.Update()

    mesh_data = PanelMesh(file_name)
    mesh_data.orig_type = 'VTK Polydata (.vtp)'
    readerOut = reader.GetOutput()
    mesh_data.num_faces =    int(readerOut.GetNumberOfCells())
    mesh_data.num_points =   int(readerOut.GetNumberOfPoints())

    for i in xrange(mesh_data.num_points):
        mesh_data.points.append(readerOut.GetPoint(i))
    mesh_data.points = np.array(mesh_data.points)

    for i in xrange(mesh_data.num_faces):
        c = readerOut.GetCell(i)
        numCellPoints = int(c.GetNumberOfPoints())
        idsTemp = []
        for i in xrange(numCellPoints):
            idsTemp.append(int(c.GetPointId(i)))
        mesh_data.faces.append(np.array(idsTemp))


    print 'Successfully read VTP mesh ' + str(file_name)

    return mesh_data


def _read_nemoh(file_name):
    '''
    Function to read nemoh meshes

    Inputs:
        file_name: name of the mesh file to be read
    Outputs:
        mesh_data: panelMesh object containing the mesh data
    '''

    with open(file_name,'r') as fid:

        lines = fid.readlines()
    temp = np.array([np.array(str(lines[i]).split()).astype(float) for i in range(1,np.size(lines))])
    count = 0

    mesh_data = PanelMesh(file_name)
    mesh_data.orig_type = 'NEMOH (.dat)'

    while temp[count,0] != 0.:

        mesh_data.points.append(temp[count,1:])
        count += 1
    count += 1
    while sum(temp[count,:]) != 0.:

        mesh_data.faces.append(temp[count,:])
        count += 1
    mesh_data.points = np.array(mesh_data.points)
    mesh_data.faces = np.array(mesh_data.faces)-1
    mesh_data.num_points = np.shape(mesh_data.points)[0]
    mesh_data.num_faces = np.shape(mesh_data.faces)[0]

    print 'Successfully read NEMOH mesh ' + str(file_name)

    return mesh_data

def _mk_vtk_id_list(it):
    '''
    Function to make vtk id list object

    Inputs:
        it: list of nodes for the face
    Outputs:
        vil: vtk id list
    '''
    vil = vtk.vtkIdList()

    for i in it:

        vil.InsertNextId(int(i))

    return vil


#def removeSurfacePanels(self):
#    tempFaces = []
#    count = 0
#
#    for i in xrange(self.num_faces):
#        deleteFace = 0
#        p0 = self.faces[i][0]
#        p1 = self.faces[i][1]
#        p2 = self.faces[i][2]
#        p3 = self.faces[i][3]
#        z0 = float(self.cords[int(p0)][2])
#        z1 = float(self.cords[int(p1)][2])
#        z2 = float(self.cords[int(p2)][2])
#        z3 = float(self.cords[int(p3)][2])
#
#        if z0 == 0.:
#            deleteFace += 1
#        if z1 == 0.:
#            deleteFace += 1
#        if z2 == 0.:
#            deleteFace += 1
#        if z3 == 0.:
#            deleteFace += 1
#        if deleteFace != 4:
#            tempFaces.append(self.faces[i])
#            count  += 1
#    print 'removed ' + str(count) + ' surface faces'
#    self.faces = []
#    self.faces = tempFaces
#    self.num_faces = np.shape(self.faces)[0]
#    os.remove(self.files['MeshInputVtp'])
#    self.write_nemohMeshInputVtp()

def read(file_name):
    file_name = os.path.abspath(file_name)
    (f_name,f_ext) = os.path.splitext(file_name)

    if f_ext == '.GDF' or f_ext == '.gdf':

        mesh_data = _read_gdf(file_name)


    elif f_ext == '.stl':

        mesh_data = _read_stl(file_name)


    elif f_ext == '.vtp':

        mesh_data = _read_vtp(file_name)


    elif f_ext == '.dat':

        mesh_data = _read_nemoh(file_name)

    else:
        raise Exception(f_ext + ' is an unsupported file mesh file type')

    mesh_data._create_vtp_mesh()

    return mesh_data
