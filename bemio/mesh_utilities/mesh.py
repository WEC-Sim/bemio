"""
Copyright 2014 the National Renewable Energy Laboratory

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

"""
This module serves the purpose of reading and writing  surface meshes from the
following mesh formats:
    * WAMTI
    * Nemoh
    * VTK Polydata (*.vtp)
        
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
    ''' Class to handel surface meshes. The class currently handels quad or tri
    meshes
    '''
    def __init__(self,file_name):

        self.file_name = file_name
        self.faces = []
        self.points = []
        self.num_faces = None
        self.num_points = None
        self.orig_type = None
        
        if os.path.isfile(file_name) is False:

            raise Exception('The file ' + file_name + ' does not exist')


    def __repr__(self):

        out_string = 'File name: ' + str(self.file_name) + \
        '\nNumber of points: ' + str(self.num_points) + \
        '\nNumber of faces: ' + str(self.num_faces) + \
        '\nOriginal mesh type: ' + str(self.orig_type) + \
        '\nObject type: bemio.mesh_utilities.mesh.PanelMesh'

        return out_string

        

    def write_vtp(self,outputFile=None):
        '''
        Function to write VTK Poly data file with the mesh data
        
        Inputs:
            outputFile: string specifying the file name for the output vtp file
        Outputs:
            None
        Function action:
            A vtp file is written with the name of the outputFile
        '''

        mesh    = vtk.vtkPolyData()
        points  = vtk.vtkPoints()
        polys   = vtk.vtkCellArray()
        
        for i in range(self.num_points):
            points.InsertPoint(i, self.points[i])
        for i in range(self.num_faces):
            polys.InsertNextCell( mk_vtk_id_list(self.faces[i]))
            
        mesh.SetPoints(points)
        mesh.SetPolys(polys)
        
        writer = vtk.vtkXMLPolyDataWriter()
        if outputFile is None:
            writer.SetFileName(self.file_name[:-3]+'vtp')
        else:
            writer.SetFileName(outputFile)
        writer.SetInputData(mesh)
        writer.SetDataModeToAscii()
        writer.Write()
        self.vtpOutFile = outputFile
        
    def write_nemoh(self,outputFile=None):
        '''
        Function to write a mesh file in the Nemoh format.
        
        The function currently assumes that the mesh that is that is read has
        four nodes per cell.
        
        Inputs:
            outputFile: string specifying the file name for the output file
        Outputs:
            None
        Function action:
            A Nemoh mesh file is written with the name of the outputFile
            
        '''

        if outputFile is None:
            with open(self.file_name[:-3]+'nemoMesh.dat','w') as fid:
                self._write_nemoh(fid)
        else:
            with open(outputFile,'w') as fid:
                self._write_nemoh(fid)
        self.nemohMeshOutFile = outputFile
                
    def _write_nemoh(self,fid):

        fid.write('2 0') # This should not be hardcoded
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
        
    def write_gdf(self,outputFile=None):
        '''
        Function to write a mesh file in the WAMIT format.
        
        The function currently assumes that the mesh that is that is read has
        four nodes per cell.
        
        Inputs:
            outputFile: string specifying the file name for the output file
        Outputs:
            None
        Function action:
            A WAMIT mesh file is written with the name of the outputFile
            
        '''

        if outputFile is None:
            with open(self.file_name[:-3]+'gdf','w') as fid:
                self._write_gdf(fid)
        else:
            with open(outputFile,'w') as fid:
                self._write_gdf(fid)
        self.gdfOutFile = outputFile

        print 'Wrote gdf file to ' + str(outputFile)
                
    def _write_gdf(self,fid):

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



    def paraview(self):
        '''
        Visualize the geometry in paraview

        To use this function make a symbolic link to the paraview.app folder
        on your system to ~/bin. Or alternatively change this function
        '''

        file_name = self.file_name[:-3] + 'vis-temp.vtp'
        if _platform == 'darwin':

            self.write_vtp(outputFile=file_name)
            try:

                os.system('open ~/bin/paraview ' + file_name + ' &')
            except:

                raise Exception('~/bin/paraview not found')
        else:

            print 'paraview() function only supported for osx'


def read(file_name):
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

    return mesh_data


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
    
    return mesh_data
    
def mk_vtk_id_list(it):
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
