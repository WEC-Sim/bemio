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
    
    def __init__(self,meshFileName):

        self.meshFileName = meshFileName
        self.faces = []
        self.points = []
        self.numFaces = None
        self.numPoints = None
        
        if os.path.isfile(meshFileName) is False:

            raise Exception('The file ' + meshFileName + ' does not exist')
        

    def writeVtp(self,outputFile=None):
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
        
        for i in range(self.numPoints):
            points.InsertPoint(i, self.points[i])
        for i in range(self.numFaces):
            polys.InsertNextCell( mkVtkIdList(self.faces[i]))
            
        mesh.SetPoints(points)
        mesh.SetPolys(polys)
        
        writer = vtk.vtkXMLPolyDataWriter()
        if outputFile is None:
            writer.SetFileName(self.meshFileName[:-3]+'vtp')
        else:
            writer.SetFileName(outputFile)
        writer.SetInputData(mesh)
        writer.SetDataModeToAscii()
        writer.Write()
        self.vtpOutFile = outputFile
        
    def writeNemoh(self,outputFile=None):
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
            with open(self.meshFileName[:-3]+'nemoMesh.dat','w') as fid:
                self._writeNemoh(fid)
        else:
            with open(outputFile,'w') as fid:
                self._writeNemoh(fid)
        self.nemohMeshOutFile = outputFile
                
    def _writeNemoh(self,fid):

        fid.write('2 0') # This should not be hardcoded
        fid.write('\n')
        for i in xrange(self.numPoints):
            fid.write(str(i+1) + ' ' +str(self.points[i]).replace('[','').replace(']',''))
            fid.write('\n')
        fid.write('0 0 0 0')
        fid.write('\n')
        for i in xrange(self.numFaces):
            fid.write(str(self.faces[i]+1).replace('[','').replace(']','').replace('.',''))
            fid.write('\n')
        fid.write('0 0 0 0')
        
    def writeGdf(self,outputFile=None):
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
            with open(self.meshFileName[:-3]+'gdf','w') as fid:
                self._writeGdf(fid)
        else:
            with open(outputFile,'w') as fid:
                self._writeGdf(fid)
        self.gdfOutFile = outputFile

        print 'Wrote gdf file to ' + str(outputFile)
                
    def _writeGdf(self,fid):

        fid.write('Mesh file written by meshio.py')
        fid.write('\n')
        fid.write('1 9.80665       ULEN GRAV')
        fid.write('\n')
        fid.write('0  0    ISX  ISY')
        fid.write('\n')
        fid.write(str(self.numFaces))
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

        fileName = self.meshFileName[:-3] + 'vis-temp.vtp'
        if _platform == 'darwin':

            self.writeVtp(outputFile=fileName)
            try:

                os.system('open ~/bin/paraview ' + fileName + ' &')
            except:

                raise Exception('~/bin/paraview not found')
        else:

            print 'paraview() function only supported for osx'

def readGdf(fileName):
    '''
    Function to read WAMIT GDF meshes
    
    Inputs:
        fileName: name of the mesh file to be read
    Outputs:
        meshData: panelMesh object containing the mesh data
    '''

    with open(fileName,'r') as fid:

        lines = fid.readlines()
        
    meshData = PanelMesh(fileName)
        
    meshData.gdfLines = lines
    meshData.uLen = int(lines[1].split()[0])  
    meshData.gravity = float(lines[1].split()[1])
    meshData.isx = float(lines[2].split()[0])
    meshData.isy = float(lines[2].split()[1])
    meshData.numFaces = int(lines[3].split()[0])
    meshData.numPoints = meshData.numFaces * 4
    meshData.points = np.array([temp.split() for temp in lines[4:]]).astype(np.float)

    meshData.pointsString = [str(temp).replace(","      ,'').replace('\r','') for temp in lines[4:]] # Output string for Nemoh mesh fil

    for panelNum,i in enumerate(np.arange(4,4+meshData.numPoints,4)):

        meshData.faces.append(np.array([i-4,i-3,i-2,i-1]))
        
    return meshData


def readStl(fileName):
    '''
    Function to read STL meshes
    
    Inputs:
        fileName: name of the mesh file to be read
    Outputs:
        meshData: panelMesh object containing the mesh data
    '''

    reader = vtk.vtkSTLReader()
    reader.SetFileName(fileName)
    reader.Update()
    
    meshData = PanelMesh(fileName)
    meshData.numFaces = int(reader.GetOutput().GetNumberOfCells())
    meshData.numPoints = meshData.numFaces * 3
    for i in range(meshData.numFaces):
        n = i*3
        meshData.faces.append(np.array([n,n+1,n+2,n+2]))
        meshData.points.append(np.array(vtk_to_numpy(reader.GetOutput().GetCell(i).GetPoints().GetData())))
    meshData.points = np.array(meshData.points).reshape([meshData.numFaces*3,3])
    
    return meshData
    

def readVtp(fileName):
    '''
    Function to read VTK Polydata meshes
    
    Inputs:
        fileName: name of the mesh file to be read
    Outputs:
        meshData: panelMesh object containing the mesh data
    '''

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(fileName)
    reader.Update()
    
    meshData = PanelMesh(fileName)
    readerOut = reader.GetOutput()
    meshData.numFaces =    int(readerOut.GetNumberOfCells())
    meshData.numPoints =   int(readerOut.GetNumberOfPoints())
 
    for i in xrange(meshData.numPoints):
        meshData.points.append(readerOut.GetPoint(i))
    meshData.points = np.array(meshData.points)

    for i in xrange(meshData.numFaces):
        c = readerOut.GetCell(i)
        numCellPoints = int(c.GetNumberOfPoints())
        idsTemp = []
        for i in xrange(numCellPoints):
            idsTemp.append(int(c.GetPointId(i)))
        meshData.faces.append(np.array(idsTemp))
    
    return meshData
    

def readNemoh(fileName):
    '''
    Function to read nemoh meshes
    
    Inputs:
        fileName: name of the mesh file to be read
    Outputs:
        meshData: panelMesh object containing the mesh data
    '''

    with open(fileName,'r') as fid:

        lines = fid.readlines()
    temp = np.array([np.array(str(lines[i]).split()).astype(float) for i in range(1,np.size(lines))])
    count = 0
    
    meshData = PanelMesh(fileName)
    
    while temp[count,0] != 0.:

        meshData.points.append(temp[count,1:])
        count += 1
    count += 1
    while sum(temp[count,:]) != 0.:

        meshData.faces.append(temp[count,:])
        count += 1
    meshData.points = np.array(meshData.points)
    meshData.faces = np.array(meshData.faces)-1
    meshData.numPoints = np.shape(meshData.points)[0]
    meshData.numFaces = np.shape(meshData.faces)[0]
    
    return meshData
    
def mkVtkIdList(it):
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
#    for i in xrange(self.numFaces):
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
#    self.numFaces = np.shape(self.faces)[0]
#    os.remove(self.files['MeshInputVtp'])
#    self.writeNemohMeshInputVtp()