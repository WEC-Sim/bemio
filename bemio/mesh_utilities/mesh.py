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
"""
This module serves provides utilities to work with the following mesh types:
    * WAMIT
    * Nemoh
    * VTK Polydata (*.vtp)
    * STL files

The functionality provided includes:
    * Ability to read mesh the mesh formats listed above
    * Ability to convert between mesh formats listed above
    * Utilities to calculate volume, surface area, linear sprint stiffness, and
    other related mesh parameters.

.. Note::
"""
import numpy as np

import imp

import os

from sys import platform as _platform

from copy import copy

from platform import system as _system

try:

    import vtk
    from vtk.util.numpy_support import vtk_to_numpy


except:

    print('The VTK Python module is required for a significant amnount of functionality in this module. Many functions will not be available for use.')



class VTK_Exception(Exception):
    pass

class PanelMesh(object):
    ''' Class to store mesh data. All mesh data is currently read and stored
    as quad elements. Tri elements are supported, but are stored as quad
    elements with a repeated point.

    Parameters:
        file_name : str
            Name of mesh file to read. Currently WAMIT (.gdf), Stereolithography
            (.stl), VTK PolyDATA (.vtp), and NEMOH (.dat) mesh formats are
            supported

    Attribuites:
        files : dictionary
            Dictionary containg input and output file names
        orig_type : str
            Mesh type of input file
        points : list
            List of points that define the mesh. `points[n] = [x coord, y coord, z coord]`.
        faces : list
            List of points that define connectivity for each face. `face[n] = [point 1 index
            , point 2 index, point 3 index, point 4 index]`, where 'point 1-4'
            are integres that correspond to the point index in the `points`
            attribuite.
        center_of_gravity : np.array
            Center of gravity of floating body
        center_of_buoyancy : np.array
            Center of buoyancy
        volume_vtk : np.array
            Mesh volume determined using VTK
        volume_x, y, and z : np.array
            Mesh volume determined using internal bemio calculations
        surface_area_vtk : float
            Surface area determined using VTK
        surface_area : float
            Surface area determined using internal bemio calculations
        normals : np.array
            Cell normals. This arrays is size `[faces.shape[0], 3]`.
            `normals[n] = [x, y, z]` is a vector that defines the normal vector
            for face `n`
        cell_surface_area : np.array
            Cell survace area. This array is size `[faces.shape[0], 3]`.
            `cell_surface_area[n]` is the surface aras of face[n].
        centroid : np.array
            Cell centroid. This array is size `[faces.shape[0], 3]`.
            `cell_surface_area[n]` is centroid of face[n].
        hydrostatic_stiffness : np.array
            The linear hydrostatic stiffness matrix of the mesh assuming the
            water surface is at z=0
        bounds : dictionary
            The bounds of the mesh. `bounds['min']` and `bounds['max']` are the
            minimum and maximum mesh dimensions, respectively
    '''
    def __init__(self,file_name):

        self.files = {}
        self.files['input_file'] = file_name

        self.orig_type = None
        self.points = []
        self.faces = []
        self.center_of_gravity = np.array([0., 0., 0.])
        self._center_of_buoyancy = None
        self._volume_vtk = None
        self._volume_x = None
        self._volume_y = None
        self._volume_z = None
        self._surface_area = None
        self._surface_area_vtk = None
        self._normals = None
        self._cell_surface_area = None
        self._centroid = None
        self._hydrostatic_stiffness = None
        self._bounds = None

        self.zero_tol = -1.e-3

        try:

            imp.find_module('vtk')
            self.VTK_installed = True

        except :

            self.VTK_installed = False

        if os.path.isfile(file_name) is False:

            raise Exception('The file ' + file_name + ' does not exist')

    def __repr__(self):
        out_string = 'Object type: bemio.mesh_utilities.mesh.PanelMesh' + \
        '\nFile name: ' + str(self.files['input_file']) + \
        '\nNumber of points: ' + str(self.points.shape[0]) + \
        '\nNumber of faces: ' + str(self.faces.shape[0]) + \
        '\nOriginal mesh type: ' + str(self.orig_type) + \
        '\nMesh bounds:' + \
        '\n\tMax: ' + str(self.bounds['max']) + \
        '\n\tMin: ' + str(self.bounds['min']) + \
        '\nCenter of mass: ' + str(self.center_of_gravity) + \
        '\nCenter of buoyancy: ' + str(self.center_of_buoyancy) + \
        '\nMesh volume [volume_x, volume_y, volume_z]: [' + str(self.volume_x) + ', ' + str(self.volume_y) + ', ' + str(self.volume_z) + ']' + \
        '\nMesh surface area: ' + str(self.surface_area) + \
        '\nHydrostatic stiffness: ' + \
        '\n\tC[3,3], C[3,4], C[3,5]: '  +  str(self.hydrostatic_stiffness[2,2]) + ', ' + str(self.hydrostatic_stiffness[2,3]) + ', ' + str(self.hydrostatic_stiffness[2,4]) + \
        '\n\tC[4,4], C[4,5], C[4,6]: '  +  str(self.hydrostatic_stiffness[3,3]) + ', ' + str(self.hydrostatic_stiffness[3,4]) + ', ' + str(self.hydrostatic_stiffness[3,5]) + \
        '\n\tC[5,5], C[5,6]:         '  +  str(self.hydrostatic_stiffness[4,4]) + ', ' + str(self.hydrostatic_stiffness[4,5])

        return out_string

    @ property
    def bounds(self, ):
        if self._bounds is None:
            self._bounds = {}
            self._bounds['max'] = self.points.max(axis=0)
            self._bounds['min'] = self.points.min(axis=0)

        return self._bounds

    @property
    def hydrostatic_stiffness(self, ):
        '''Getter for the `hydrostatic_stiffness` variable.

        Calculated as defined in Section 3.1 of the WAMIT v7.0 users manual.
        '''
        if self._hydrostatic_stiffness is None:

            self._hydrostatic_stiffness = np.zeros([6,6])
            for face_n,face in enumerate(self.faces):

                if self.points[face[0]][2] <= self.zero_tol or \
                self.points[face[1]][2] <= self.zero_tol or \
                self.points[face[2]][2] <= self.zero_tol or \
                self.points[face[3]][2] <= self.zero_tol:
                    self._hydrostatic_stiffness[2,2] += -self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[2,3] += -self.centroid[face_n][1] * self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[2,4] += self.centroid[face_n][0] * self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[3,3] += -self.centroid[face_n][1]**2 * self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[3,4] += self.centroid[face_n][0] * self.centroid[face_n][1] * self.normals[face_n][2] * self.cell_surface_area[face_n]
                    self._hydrostatic_stiffness[4,4] += -self.centroid[face_n][0]**2 * self.normals[face_n][2] * self.cell_surface_area[face_n]

            self._hydrostatic_stiffness[3,3] += self.volume_x * self.center_of_buoyancy[2] - self.volume_x * self.center_of_gravity[2]
            self._hydrostatic_stiffness[3,5] += -self.volume_x * self.center_of_buoyancy[0] + self.volume_x * self.center_of_gravity[0]
            self._hydrostatic_stiffness[4,4] += self.volume_x * self.center_of_buoyancy[2] - self.volume_x * self.center_of_gravity[2]
            self._hydrostatic_stiffness[4,5] += -self.volume_x * self.center_of_buoyancy[1] + self.volume_x * self.center_of_gravity[1]

            print('Calculated hydorstatic stiffness')

        return self._hydrostatic_stiffness

    @property
    def center_of_buoyancy(self, ):
        '''Getter for the `center_of_buoyancy` variable.

        Calculated as defined in Section 3.1 of the WAMIT v7.0 users manual.
        '''
        if self._center_of_buoyancy is None:
            x_b = 0.
            y_b = 0.
            z_b = 0.
            self._center_of_buoyancy = 0.

            for face_n,face in enumerate(self.faces):
                if self.points[face[0]][2] <= self.zero_tol or \
                self.points[face[1]][2] <= self.zero_tol or \
                self.points[face[2]][2] <= self.zero_tol or \
                self.points[face[3]][2] <= self.zero_tol:
                    x_b += self.normals[face_n][0]*self.centroid[face_n][0]**2*self.cell_surface_area[face_n]
                    y_b += self.normals[face_n][1]*self.centroid[face_n][1]**2*self.cell_surface_area[face_n]
                    z_b += self.normals[face_n][2]*self.centroid[face_n][2]**2*self.cell_surface_area[face_n]
            self._center_of_buoyancy = 1./(2.*self.volume_x )*np.array([x_b, y_b, z_b])

            print('Calculated the center of buoyancy')

        return self._center_of_buoyancy

    @property
    def normals(self, ):
        if self._normals is None:
            self._normals = {}
            for face_n in range(self.faces.shape[0]):
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

            print('Calculated mesh cell normals')

        return self._normals

    @property
    def cell_surface_area(self):
        '''Getter for `cell_surface_area`

        Calculated
        '''
        if self._cell_surface_area is None:
            self._cell_surface_area = {}
            for face_n in range(self.faces.shape[0]):
                a = self.points[self.faces[face_n][1]] - self.points[self.faces[face_n][0]]
                b = self.points[self.faces[face_n][2]] - self.points[self.faces[face_n][1]]
                c = self.points[self.faces[face_n][3]] - self.points[self.faces[face_n][2]]
                d = self.points[self.faces[face_n][0]] - self.points[self.faces[face_n][3]]
                self._cell_surface_area[face_n] = 1./2. * ( np.linalg.norm(np.cross(a,b)) + np.linalg.norm(np.cross(c,d)) )

            print('Calculated surface cell area')

        return self._cell_surface_area

    @property
    def surface_area(self):
        if self._surface_area is None:

            self._surface_area = sum(self.cell_surface_area.values())

            print('Calculated surface area')

        return self._surface_area

    @property
    def volume_vtk(self):
        if self.VTK_installed is False:
            raise VTK_Exception('VTK must be installed to access the volume_vtk property')

        if self._volume_vtk is None:
            tri_converter = vtk.vtkTriangleFilter()
            tri_converter.SetInputDataObject(self.vtp_mesh)
            tri_converter.Update()
            tri_mesh = tri_converter.GetOutput()
            mass_props = vtk.vtkMassProperties()
            mass_props.SetInputDataObject(tri_mesh)
            self._volume_vtk = mass_props.GetVolume()

            print('Calculated mesh volume using VTK library')

        return self._volume_vtk

    @property
    def centroid(self):

        if self._centroid is None:
            self._centroid = {}

            for face_n,face in enumerate(self.faces):
                points = [self.points[face[0]],
                          self.points[face[1]],
                          self.points[face[2]],
                          self.points[face[3]]]

                points = list(map(np.asarray, set(map(tuple, points)))) # This removes duplicate points... somehow
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

    @property
    def surface_area_vtk(self):
        if self.VTK_installed is False:
            raise VTK_Exception('VTK must be installed to access the surface_area_vtk property')
        if self._surface_area_vtk is None:
            tri_converter = vtk.vtkTriangleFilter()
            tri_converter.SetInputDataObject(self.vtp_mesh)
            tri_converter.Update()
            tri_mesh = tri_converter.GetOutput()
            mass_props = vtk.vtkMassProperties()
            mass_props.SetInputDataObject(tri_mesh)
            self._surface_area_vtk = mass_props.GetSurfaceArea()

            print('Calculated mesh surface area using VTK Python bindings')

        return self._surface_area_vtk

    def write(self,mesh_format='VTP'):
        '''Function to write NEMOH, WAMIT, or VTK PolyData formats.

        Parameters:
            mesh_format : string {'VTP', 'WAMIT', 'NEMOH'}
                Variable that specifies the mesh format to write.

        Examples:
            This example assumes that a mesh has been read by bemio and mesh
            data is contained in a `PanelMesh` object called `mesh`

            Here is how a WAMIT mesh would be written
            >>> mesh.write(mesh_format='WAMTI')
        '''
        if mesh_format == 'VTK' or mesh_format == 'VTP':
            self._write_vtp()

        if mesh_format == 'WAMIT' or mesh_format == 'GDF':
            self._write_gdf()

        if mesh_format == 'NEMOH':
            self._write_nemoh()

    def open(self):
        '''Function to open a VTK PolyData object in the default viewer of your
        operating system.

        .. Note::
            This function is only available for OSX and Linux systems and
            and requires you have a program installed that has the ability to
            open VTK PolyData (.vtp) files.

        Example:
            This example assumes that a mesh has been read by bemio and mesh
            data is contained in a `PanelMesh` object called `mesh`

            >>> mesh.open()
        '''
        self.write(mesh_format='VTP')
        if _system() == 'Darwin':
            os.system('open ' + self.files['vtp'])

        elif _system() == 'Linux':
            os.system('xdg ' + self.files['vtp'])

        else:
            raise Exception('The open function is only supported for OSX')

    def calculate_center_of_gravity_vtk(self, ):
        '''Function to calculate the center of gravity

        .. Note::
            The VTK Pytnon bindings must be installed to use this function

        Examples:
            This example assumes that a mesh has been read by bemio and mesh
            data is contained in a `PanelMesh` object called `mesh`

            >>> mesh.calculate_center_of_gravity_vtk()
        '''
        if self.VTK_installed is False:
            raise VTK_Exception('VTK must be installed to access the calculate_center_of_gravity_vtk function')

        com = vtk.vtkCenterOfMass()
        if vtk.VTK_MAJOR_VERSION >= 6:
            com.SetInputData(self.vtp_mesh)
        else:
            com.SetInput(self.vtp_mesh)
        com.Update()
        self.center_of_gravity = com.GetCenter()

        print('Calculated center of gravity assuming uniform material density')

    def view(self, color=[0.5,1,0.5], opacity=1.0, save_png=False, camera_pos=[50,50,50], interact=True):
        '''Function to view the mesh using the VTK library

        Parameters:
            color : list, optional
                VTK color specification for the mesh
            opackty : float, optional
                VTK opacity for the mesh. Must be between 0. and 1.
            save_png : bool
                Boolean operater that determines if a .png image of the mesh is
                saved.
            interact : bool, optional
                Boolean operater that determines if the user can interact with
                the geometry (e.g. zoom and rotate) after it is displayed
            camera_pos : list, optional

        Examples:
            This example assumes that a mesh has been read by bemio and mesh
            data is contained in a `PanelMesh` object called `mesh`

            >>> mesh.view()

            The mesh view window must be closed in order to return command to
            the Python shell
        '''
        if self.VTK_installed is False:
            raise VTK_Exception('VTK must be installed to use the view function')

        # Create a mapper and load VTP data into the mapper
        mapper=vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION >= 6:
            mapper.SetInputData(self.vtp_mesh)
        else:
            mapper.SetInput(self.vtp_mesh)

        # Create an actor that contains the data in the mapper
        actor=vtk.vtkActor()
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetOpacity(opacity)
        actor.SetMapper(mapper)
        actor.GetProperty().EdgeVisibilityOn()

        # Camera
        camera = vtk.vtkCamera();
        camera.SetPosition(camera_pos)
        camera.SetFocalPoint(0, 0, 0)

        # Add axes
        axes = vtk.vtkAxesActor()

        # Render the data
        ren = vtk.vtkRenderer()
        ren.AddActor(actor)
        ren.AddActor(axes)
        ren.SetActiveCamera(camera)

        # Create a render window
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        renWin.SetSize(800, 800)

        # Start the visiuilization
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        ren.SetBackground(0,0,0)
        renWin.Render()


        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()

        if save_png is True:
            w2if = vtk.vtkWindowToImageFilter()
            w2if.SetInput(renWin)
            w2if.Update()

            writer = vtk.vtkPNGWriter()
            writer.SetFileName(self.files['png'])
            writer.SetInputDataObject(w2if.GetOutput())
            writer.Write()

            print('Wrote mesh image to: ' + self.files['png'])

        if interact is True:
            iren.Start()

    def view_points_and_vectors(self, point_show=True, point_size=7.5, point_color=[1., 0., 0.], point_opacity=1., point_colorByScalar=False, point_scalar=[], normal_show=True, normal_length=0.3, normal_tip_radius=0.1, normal_scaleByPointScalar=False, normal_scale=1, normal_color=[0., 0., 1.], normal_opacity=1., normal_reverse=False, normal_arrowAtPoint=False, normal_colorByPointScalar=False, mesh_show=True, mesh_color=[1., 1., 1.], mesh_opacity=1, mesh_colorByScalar=False, mesh_scalar=[], camera_pos=[50,50,50], save_png=False, interact=True):
        '''Function to view the mesh using the VTK library. Allows for 
        visualization of mesh, points (centroids), and normal vectors, as well
        as scalar fields. prior to scaling

        Parameters:
            point_show : bool
                Boolean operator that determines if the points are shown.
            point_size : float
                Size of points
            point_color : list, float
                VTK color specification for points [r,g,b].
            point_opacity : float
                VTK opacity for the points. Must be between 0. and 1.
            point_colorByScalar : bool
                Boolean operator that determines if the points are colored by a
                scalar field 'point_scalar' rather than by constant 
                'point_color'
            point_scalar : np.array
                scalar value at each point.
            normal_show : bool
                Boolean operator that determines if the normals are shown.
            normal_length : float
                length of the arrows (glyphs) prior to scaling
            normal_tip_radius: float 
                Tip radius of the arrows (glyphs) prior to scaling
            normal_scaleByPointScalar : bool
                Boolean operator that determines if the normal vectors are
                scaled by the scalar field 'point_scalar'
            normal_scale : float
                Value to scale the size of the normal glyphs by
            normal_color : list, float
                VTK color specification for the normals [r,g,b].
            normal_opacity : float
                VTK opacity for the normals. Must be between 0. and 1.
            normal_reverse : bool
                Boolean operator that determines if the direction of the normal
                vectors is reversed 
            normal_arrowAtPoint: bool
                Boolean operator that determines if the arrow is at the point 
                or at the opposite edge of the glyph
            normal_colorByPointScalar : bool
                Boolean operator that determines if the normals are colored by a
                scalar field 'point_scalar' rather than by constant 
                'normal_color'
            mesh_show : bool
                Boolean operator that determines if the points are shown.
            mesh_color : list, float
                VTK color specification for the mesh [r,g,b].
            mesh_opacity : float
                VTK opacity for the mesh. Must be between 0. and 1.
            mesh_colorByScalar : bool
                Boolean operator that determines if the mesh is colored by a
                scalar field 'mesh_scalar' rather than by a constant 
                'mesh_color'
            mesh_scalar : np.array
                scalar value at each mesh cell.
            camera_pos : list, floats
                Position of camera [X Y Z].
            save_png : bool
                Boolean operater that determines if a .png image of the mesh is
                saved.
            interact : bool, optional
                Boolean operater that determines if the user can interact with
                the geometry (e.g. zoom and rotate) after it is displayed

        Examples:
            This example assumes that a mesh has been read by bemio and mesh
            data is contained in a `PanelMesh` object called `mesh`

            >>> mesh.view_points_and_vectors()

            The mesh view window must be closed in order to return command to
            the Python shell
        '''
        if self.VTK_installed is False:
            raise VTK_Exception('VTK must be installed to use the view function')
        # Centroid
        centroid_points = vtk.vtkPoints()
        vertices = vtk.vtkCellArray()
        for ip in range(len(self.centroid)):
            id = centroid_points.InsertNextPoint(list(self.centroid[ip]))
            vertices.InsertNextCell(1)
            vertices.InsertCellPoint(id)
        centroid = vtk.vtkPolyData()
        centroid.SetPoints(centroid_points)
        centroid.SetVerts(vertices)
        n=np.array(list(self.normals.values()))
        nv = vtk.util.numpy_support.numpy_to_vtk(n)
        pv =centroid.GetPointData()
        _ = pv.SetNormals(nv)
        if point_colorByScalar or normal_colorByPointScalar:
            point_scalar = point_scalar.astype('double')
            vpoint_scalar =  vtk.util.numpy_support.numpy_to_vtk(np.ascontiguousarray(point_scalar))
            _ = pv.SetScalars(vpoint_scalar)
            point_scalarRange = centroid.GetScalarRange()
        centroid.Modified()
        centroidMapper=vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION >= 6:
            centroidMapper.SetInputData(centroid)
        else:
            centroidMapper.SetInput(centroid)
        if point_colorByScalar:
            centroidMapper.SetScalarModeToUsePointData
            centroidMapper.SetColorModeToMapScalars()
            centroidMapper.SetScalarRange(point_scalarRange)
        else:
            centroidMapper.SetColorModeToDefault()
            centroidMapper.SetScalarVisibility(0)
        centroidActor = vtk.vtkActor()
        centroidActor.SetMapper(centroidMapper)
        centroidActor.GetProperty().SetPointSize(point_size)
        centroidActor.GetProperty().SetOpacity(point_opacity)
        if not point_colorByScalar:
            centroidActor.GetProperty().SetColor(point_color)
        # Normals
        def MakeGlyphs(src, normal_reverse, normal_arrowAtPoint):
            reverse = vtk.vtkReverseSense()
            maskPts = vtk.vtkMaskPoints()
            maskPts.SetOnRatio(1)
            if normal_reverse:
                reverse.SetInputData(src)
                reverse.ReverseCellsOn()
                reverse.ReverseNormalsOn()
                maskPts.SetInputConnection(reverse.GetOutputPort())
            else:
                maskPts.SetInputData(src)
            arrow = vtk.vtkArrowSource()
            if normal_arrowAtPoint:
                arrow.SetInvert(1)
            else: 
                arrow.SetInvert(0)
            arrow.SetTipResolution(16)
            arrow.SetTipLength(normal_length)
            arrow.SetTipRadius(normal_tip_radius)
            glyph = vtk.vtkGlyph3D()
            glyph.SetSourceConnection(arrow.GetOutputPort())
            glyph.SetInputConnection(maskPts.GetOutputPort())
            glyph.SetVectorModeToUseNormal()
            glyph.SetScaleFactor(normal_scale)
            if normal_scaleByPointScalar:
                glyph.SetScaleModeToScaleByScalar()
            else:
                glyph.SetScaleModeToScaleByVector()
            glyph.SetColorModeToColorByScalar()
            glyph.OrientOn()
            return glyph
        glyph = MakeGlyphs(centroid,normal_reverse,normal_arrowAtPoint)
        glyphMapper = vtk.vtkPolyDataMapper()
        glyphMapper.SetInputConnection(glyph.GetOutputPort())
        glyphActor = vtk.vtkActor()
        glyphActor.SetMapper(glyphMapper)
        glyphActor.GetProperty().SetOpacity(normal_opacity)
        if not normal_colorByPointScalar:
            glyphActor.GetProperty().SetColor(normal_color)
        if normal_colorByPointScalar:
            glyphMapper.SetScalarModeToUsePointData()
            glyphMapper.SetScalarVisibility(1)
            glyphMapper.SetColorModeToMapScalars()
            normal_scalarRange = centroid.GetScalarRange()
            glyphMapper.SetScalarRange(normal_scalarRange)
        else:
            glyphMapper.SetScalarVisibility(0)
            glyphMapper.SetColorModeToDefault()
        # Mesh
        if mesh_colorByScalar:
            cv = self.vtp_mesh.GetCellData()
            mesh_scalar = mesh_scalar.astype('double')
            vmesh_scalar =  vtk.util.numpy_support.numpy_to_vtk(np.ascontiguousarray(mesh_scalar))
            _ = cv.SetScalars(vmesh_scalar)
        meshMapper=vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION >= 6:
            meshMapper.SetInputData(self.vtp_mesh)
        else:
            meshMapper.SetInput(self.vtp_mesh)
        if mesh_colorByScalar:
            meshMapper.SetScalarModeToUseCellData
            meshMapper.SetColorModeToMapScalars()
            mesh_scalarRange = self.vtp_mesh.GetScalarRange()
            meshMapper.SetScalarRange(mesh_scalarRange)
        else:
            meshMapper.SetScalarVisibility(0)
        meshActor=vtk.vtkActor()
        meshActor.SetMapper(meshMapper)
        meshActor.GetProperty().EdgeVisibilityOn()
        meshActor.GetProperty().SetOpacity(mesh_opacity)
        if not mesh_colorByScalar:
            meshActor.GetProperty().SetColor(mesh_color)
        # bars
        if point_colorByScalar or normal_colorByPointScalar:
            lut = centroidMapper.GetLookupTable()
            lut.SetNumberOfTableValues(256)
            colorTransferFunction = vtk.vtkColorTransferFunction()
            colorTransferFunction.SetColorSpaceToDiverging()
            colorTransferFunction.AddRGBPoint(0,0.231373,0.298039,0.752941)
            colorTransferFunction.AddRGBPoint(0.5,0.865003,0.865003,0.865003)
            colorTransferFunction.AddRGBPoint(1.0,0.705882,0.0156863,0.14902)
            for ii,ss in enumerate([float(xx)/float(256) for xx in range(256)]):
                cc = colorTransferFunction.GetColor(ss)
                lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
            lut.Build()
            centroidMapper.SetLookupTable(lut)
            glyphMapper.SetLookupTable(lut)
            point_bar = vtk.vtkScalarBarActor()
            point_bar.SetLookupTable(lut)
            point_bar.SetTitle('Point Scalar')
        if mesh_colorByScalar:
            mesh_lut = meshMapper.GetLookupTable()
            mesh_lut.SetNumberOfTableValues(256)
            mesh_colorTransferFunction = vtk.vtkColorTransferFunction()
            mesh_colorTransferFunction.SetColorSpaceToDiverging()
            mesh_colorTransferFunction.AddRGBPoint(0,0.231373,0.298039,0.752941)
            mesh_colorTransferFunction.AddRGBPoint(0.5,0.865003,0.865003,0.865003)
            mesh_colorTransferFunction.AddRGBPoint(1.0,0.705882,0.0156863,0.14902)
            for ii,ss in enumerate([float(xx)/float(256) for xx in range(256)]):
                mesh_cc = mesh_colorTransferFunction.GetColor(ss)
                mesh_lut.SetTableValue(ii,mesh_cc[0],mesh_cc[1],mesh_cc[2],1.0)
            mesh_lut.Build()
            meshMapper.SetLookupTable(mesh_lut)
            mesh_bar = vtk.vtkScalarBarActor()
            mesh_bar.SetLookupTable(mesh_lut)
            mesh_bar.SetTitle('Mesh Scalar')
        # Camera
        camera = vtk.vtkCamera();
        camera.SetPosition(camera_pos)
        camera.SetFocalPoint(0, 0, 0)
        # Add axes
        axes = vtk.vtkAxesActor()
        # Render the data
        ren = vtk.vtkRenderer()
        if mesh_show:
            ren.AddActor(meshActor)
        if point_show:
            ren.AddActor(centroidActor)
        if normal_show:
            ren.AddActor(glyphActor)
        ren.AddActor(axes)
        ren.SetActiveCamera(camera)
        ren.SetBackground(0,0,0)
        # Create a render window
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        renWin.SetSize(800, 800)
        # Create interactive renderer
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        # Start the bar widgets
        if point_colorByScalar or normal_colorByPointScalar:
            point_bar_widget = vtk.vtkScalarBarWidget()
            point_bar_widget.SetInteractor(iren)
            point_bar_widget.SetScalarBarActor(point_bar)
            point_bar_widget.On()
        if mesh_colorByScalar:
            mesh_bar_widget = vtk.vtkScalarBarWidget()
            mesh_bar_widget.SetInteractor(iren)
            mesh_bar_widget.SetScalarBarActor(mesh_bar)
            mesh_bar_widget.On()
        # Render
        renWin.Render()
        vtk.vtkPolyDataMapper().SetResolveCoincidentTopologyToPolygonOffset()
        # Save image
        if save_png is True:
            w2if = vtk.vtkWindowToImageFilter()
            w2if.SetInput(renWin)
            w2if.Update()
            writer = vtk.vtkPNGWriter()
            writer.SetFileName(self.files['png'])
            writer.SetInputDataObject(w2if.GetOutput())
            writer.Write()
            print('Wrote mesh image to: ' + self.files['png'])
        # Set interaction mode
        if interact is True:
            iren.Start()

    def cut(self,plane=2,value=0.0,direction=1): #NOT IMPLEMENTED 
        '''This function is not currently working 100%
        '''
        raise NotImplementedError()
        #self.collapse(plane,value,direction)
        #
        #tempFaces = []
        #count = 0
        #
        #for i in xrange(self.faces.shape[0]):
        #
        #   delete_face = 0
        #
        #   for j in xrange(4):
        #
        #       p = self.faces[i][j]
        #       z = float(self.cords[int(p)][2])
        #
        #       if z == 0.:
        #           delete_face += 1
        #
        #   if delete_face != 4:
        #       tempFaces.append(self.faces[i])
        #       count  += 1
        #
        #print 'removed ' + str(count) + ' surface faces'
        #self.faces = tempFaces
        #self.faces.shape[0] = self.faces.shape[0]

    def scale(self, scale_vect):
        '''Function used to scale mesh objects in the x, y, and z directions.

        Parameters:
            scale_vect : list
                A list that contains the x, y, and z scale factors for the mesh

        Examples:
            This example assumes that a mesh has been read by bemio and mesh
            data is contained in a `PanelMesh` object called `mesh`

            Here is how to scale a mesh by a factor of 2 in the x direction and
            .5 in the y direction:

            >>> mesh.scale(scale_vect=[2, 0.5, 1])
        '''
        scale_vect = np.array(scale_vect)
        if scale_vect.size != 3:
            raise Exception('The scale_vect input must be a length 3 vector')
        self.points = self.points*scale_vect
        self.scale_vect = scale_vect
        self._create_vtp_mesh()
        print('Scaled mesh by: ' + str(scale_vect))

    def translate(self,translation_vect,translate_cog=True):
        '''Function used to translate mesh obvjects in the x, y, and z directions.

        Parameters:
            translation_vect : list
                A list that contains the desired x, y, and z translation for the
                mesh

        Examples:
            This example assumes that a mesh has been read by bemio and mesh
            data is contained in a `PanelMesh` object called `mesh`

            Here is how to translate a mesh by 2 in the x direction and
            .5 in the y direction:

            >>> mesh.translate(scale_vect=[2, 0.5, 0])
        '''
        translation_vect = np.array(translation_vect)
        if translation_vect.size != 3:
            raise Exception('The translation_vect input must be a length 3 vector')
        self.points += translation_vect
        self.translation_vect = translation_vect

        if translate_cog is True:
            self.center_of_gravity += translation_vect

        print('Translated mesh by: ' + str(translation_vect) + '\nCenter of gravity is: ' + str(self.center_of_gravity))

    def xzmirror(self,):
        '''Function used to mirror the mesh object about the xz plane.
        '''
        npoint = np.shape(self.points)[0]
        # create symmetry points
        psym = np.copy(self.points)
        psym[:,1] = psym[:,1]*-1.
        mask = np.where(psym[:,1]==0.)
        psym = np.delete(psym,mask,0)
        # create symmetry point map
        ptmap = np.zeros([npoint,2])
        ptmap[:,0] = np.arange(0,npoint)
        count = 0
        for ipoint in range(npoint):
            if self.points[ipoint][1] == 0.:
                ptmap[ipoint,1] = ptmap[ipoint,0]
            else:
                ptmap[ipoint,1] = npoint + count
                count += 1
        # append points
        self.points = np.concatenate(([np.array(self.points)],[psym]),axis=1)[0]
        # create symmetry faces
        nface = np.shape(self.faces)[0]
        fsym = [[0,0,0,0]]
        fsym_row = [0,0,0,0]
        for iface in range(nface):
            ftmp = self.faces[iface]
            fsym_row[0] = int(ptmap[ np.where(ptmap[:,0]==ftmp[3])[0] , 1])
            fsym_row[1] = int(ptmap[ np.where(ptmap[:,0]==ftmp[2])[0] , 1])
            fsym_row[2] = int(ptmap[ np.where(ptmap[:,0]==ftmp[1])[0] , 1])
            fsym_row[3] = int(ptmap[ np.where(ptmap[:,0]==ftmp[0])[0] , 1])
            fsym.append(list(fsym_row))
        fsym = np.array(fsym)[1:,:]
        self.faces = np.concatenate(([self.faces],[fsym]),axis=1)[0]
        self.sym = 0
        self.num_points = len(self.points)
        self.num_faces = len(self.faces)
        self._create_vtp_mesh()

    def remove_duplicate_points(self, ):
        '''Function to remove duplicate points.
        '''
        def unique_rows(a):
            a = np.ascontiguousarray(a)
            unique_a,inverse_a = np.unique(a.view([('', a.dtype)]*a.shape[1]), return_inverse=True)
            unique_a =  unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
            return unique_a,inverse_a
        unique_points,inverse_points = unique_rows(self.points)
        map = np.zeros([len(inverse_points),2])
        map[:,0] = np.array(list(range(len(self.points))))
        map[:,1] = inverse_points
        faces_org = self.faces.flatten()
        faces_new = np.zeros(np.shape(faces_org))
        for ip in range(len(faces_new)):
            faces_new[ip]=map[faces_org[ip],1]
        faces_new = faces_new.reshape((len(faces_new)/4,4))
        self.points = unique_points
        self.faces = faces_new
        self.num_points = len(self.points)
        self._create_vtp_mesh()

    def _create_vtp_mesh(self):
        '''Internal function to creat a VTP mesh from the imported mesh data
        '''
        if self.VTK_installed is True:

            self.vtp_mesh    = vtk.vtkPolyData()
            points  = vtk.vtkPoints()
            polys   = vtk.vtkCellArray()

            for i in range(np.array(self.points).shape[0]):

                points.InsertPoint(i, self.points[i])


            for i in range(np.array(self.faces).shape[0]):

                polys.InsertNextCell(_mk_vtk_id_list(self.faces[i]))

            self.vtp_mesh.SetPoints(points)
            self.vtp_mesh.SetPolys(polys)

    def _collapse(self,plane=2,value=0.0,direction=1): #NOT IMPLEMENTED
        #This function is not yet working 100%
        raise NotImplementedError()
        #'''Collapse points
        #'''
        #for face,face_n in xrange(self.faces.shape[0]):
        #
        #    for j in xrange(self.faces[i].size):
        #
        #        p = int(self.faces[i][j])
        #
        #        if self.points[p][plane] > value*direction:
        #
        #            self.points[p][plane] = value

    def _write_vtp(self):
        '''Internal function to write VTK PolyData mesh files
        '''
        if self.VTK_installed is False:

            raise VTK_Exception('VTK must be installed write VTP/VTK meshes, please select a different output mesh_format')

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(self.files['vtp'])
        if vtk.VTK_MAJOR_VERSION >= 6:
            writer.SetInputData(self.vtp_mesh)
        else:
            writer.SetInput(self.vtp_mesh)
        writer.SetDataModeToAscii()
        writer.Write()

        print('Wrote VTK PolyData mesh to: ' + str(self.files['vtp']))

    def _write_nemoh(self):
        '''Internal function to write NEMOH mesh files
        '''
        with open(self.files['nemoh'],'w') as fid:
            fid.write('2 0') # This should not be hard coded
            fid.write('\n')
            for i in range(self.points.shape[0]):
                fid.write(str(i+1) + ' ' +str(self.points[i]).replace('[','').replace(']',''))
                fid.write('\n')
            fid.write('0 0 0 0')
            fid.write('\n')
            for i in range(self.faces.shape[0]):
                fid.write(str(self.faces[i]+1).replace('[','').replace(']','').replace('.',''))
                fid.write('\n')
            fid.write('0 0 0 0')

        print('Wrote NEMOH mesh to: ' + str(self.files['nemoh']))

    def _write_gdf(self):
        '''Internal function to write WAMIT mesh files
        '''
        with open(self.files['wamit'],'w') as fid:
            fid.write('Mesh file written by meshio.py')
            fid.write('\n')
            fid.write('1 9.80665       ULEN GRAV')
            fid.write('\n')
            fid.write('0  0    ISX  ISY')
            fid.write('\n')
            fid.write(str(self.faces.shape[0]))
            fid.write('\n')
            for i,face in enumerate(self.faces):
                if np.size(face) is 4: # if the mesh element is a quad
                    for j,pointKey in enumerate(face):
                        fid.write(str(self.points[pointKey]).replace(',','').replace('[','').replace(']','') + '\n')
                if np.size(face) is 3: # if the mesh element is a tri
                    faceMod = np.append(face,face[-1])
                    for j,pointKey in enumerate(faceMod):
                        fid.write(str(self.points[pointKey]).replace(',','').replace('[','').replace(']','') + '\n')

        print('Wrote WAMIT mesh to: ' + str(self.files['wamit']))

    def _calc_component_vol(self, ):
        '''Internal function to calculate mesh volume using the methods
        described in Section 3.1 the WAMIT v7.0 users manual.
        '''
        self._volume_x = 0.
        self._volume_y = 0.
        self._volume_z = 0.
        volume = 0.

        for face_n in range(self.faces.shape[0]):
            volume += self.normals[face_n]*self.centroid[face_n]*self.cell_surface_area[face_n]

        self._volume_x = volume[0]
        self._volume_y = volume[1]
        self._volume_z = volume[2]

        print('Calculated x y and z mesh volumes')

def _read_gdf(file_name):
    '''Internal function to read gdf wamit meshes
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

    mesh_data.faces = np.array(mesh_data.faces)

    return mesh_data

def _read_stl(file_name):
    '''Internal function to read stl mesh files
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
    '''Internal function to read vtp mesh files
    '''

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_name)
    reader.Update()

    mesh_data = PanelMesh(file_name)
    mesh_data.orig_type = 'VTK Polydata (.vtp)'
    readerOut = reader.GetOutput()
    mesh_data.num_faces =    int(readerOut.GetNumberOfCells())
    mesh_data.num_points =   int(readerOut.GetNumberOfPoints())

    for i in range(mesh_data.num_points):
        mesh_data.points.append(readerOut.GetPoint(i))
    mesh_data.points = np.array(mesh_data.points)

    for i in range(mesh_data.num_faces):
        c = readerOut.GetCell(i)
        numCellPoints = int(c.GetNumberOfPoints())
        idsTemp = []
        for i in range(numCellPoints):
            idsTemp.append(int(c.GetPointId(i)))
        mesh_data.faces.append(np.array(idsTemp))
    mesh_data.faces = np.array(mesh_data.faces)


    return mesh_data

def _read_nemoh(file_name):
    '''Internal function to read nemoh mesh
    '''

    with open(file_name,'r') as fid:

        lines = fid.readlines()
    sym = np.array(str(lines[0]).split()).astype(int)[1]
    temp = np.array([np.array(str(lines[i]).split()).astype(float) for i in range(1,np.size(lines))])
    count = 0

    mesh_data = PanelMesh(file_name)
    mesh_data.orig_type = 'NEMOH (.dat)'
    mesh_data.sym = sym

    while temp[count,0] != 0.:

        mesh_data.points.append(temp[count,1:])
        count += 1
    count += 1
    while sum(temp[count,:]) != 0.:

        mesh_data.faces.append(temp[count,:]-1)
        count += 1
    mesh_data.points = np.array(mesh_data.points)
    mesh_data.faces = np.array(mesh_data.faces)
    mesh_data.num_points = np.shape(mesh_data.points)[0]
    mesh_data.num_faces = np.shape(mesh_data.faces)[0]

    return mesh_data

def read(file_name):
    '''Function to read surface mesh files. Currently VTK PolyData (.vtk),
    WAMIT (.gdf), NEMOH (.dat), and Stereolithography (.stl) mesh formats are
    supported

    Parameters:
        file_name : str
            Name of the mesh file

    Returns:
        mesh_data : PanelMesh
            A PanelMesh object that contains the mesh data

    Exmaples:
        This example assumes that a VTK PlolyData mesh named mesh.vtp exists
        in the current working directory

        >>> mesh = read('mesh.vtp')

        The mesh can then be converted to another format using the `write`
        function. In this case a wamit mesh is created.

        >>> mesh.write(mesh_format='WAMIT')

        If the VTK python bindings are installed the mesh can be viewed using
        the following command:

        >>> mesh.view()

        If you would are using OSX or Linux and have Paraview installed you can
        view the file using the follwing command:

        >>> mesh.open()
    '''
    print('Reading mesh file: ' + str(file_name))

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


    mesh_data.files_base = os.path.splitext(file_name)[0] + '_bemio_output'
    mesh_data.files['vtp'] = mesh_data.files_base + '.vtp'
    mesh_data.files['wamit'] = mesh_data.files_base + '.gdf'
    mesh_data.files['nemoh'] = mesh_data.files_base + '.dat'
    mesh_data.files['png'] = os.path.splitext(file_name)[0] + '.png'

    if mesh_data.VTK_installed is True:

        mesh_data._create_vtp_mesh()

    print('Successfully read mesh file: ' + str(file_name))

    return mesh_data

def _mk_vtk_id_list(it):
    '''
    Internal function to make vtk id list object

    Parameters:
        it : list
            List of nodes that define a face
    Returns:
        vil: vtkIdList
            A vtkIdList object
    '''
    vil = vtk.vtkIdList()
    for i in it:
        vil.InsertNextId(int(i))

    return vil

def collapse_to_plane(mesh_obj, plane_ind=2, plane_loc=-1e-5, cut_dir=1.): #NOT IMPLEMENTED
    '''Function to collapse points to a given plane

    .. Note::
        This function is not yet implemented
    '''
    raise NotImplementedError()

def cut_mesh(mesh_obj, plane_ind=2, plane_loc=-1e-5, cut_dir=1.):
    '''Function to remove cells on one side of plane

    .. Note::
        This function is still early in the stages of development and needs
        to be improved and made more robust.

    Parameters:
        mesh_obj : PanelMesh
            Mesh object to cut
        plane_ind : int, optional
            Index of plane along which to cut the mesh, 0 == x, 1 == y, 2 == z
        plane_loc : float
            Location of the mesh cut
        cut_dir : int, {1, -1}
            Direction for the mesh cut

    Returns:
        cut_mesh : PanelMesh
            Panel mesh object that has been cut as specified

    Examples:
        None avaiable to data
    '''
    cut_mesh = copy(mesh_obj)
    tempFaces = []
    cut_mesh.removed_faces = []
    cut_mesh.removed_points = []

    for face_n,face in enumerate(cut_mesh.faces):
        if cut_mesh.points[face[0]][plane_ind] <= plane_loc*cut_dir or \
        cut_mesh.points[face[1]][plane_ind] <= plane_loc*cut_dir or \
        cut_mesh.points[face[2]][plane_ind] <= plane_loc*cut_dir or \
        cut_mesh.points[face[3]][plane_ind] <= plane_loc*cut_dir:

            tempFaces.append(face)

        else:
            cut_mesh.removed_faces.append(face)
            cut_mesh.removed_points.append(cut_mesh.points[face[0]])
            cut_mesh.removed_points.append(cut_mesh.points[face[1]])
            cut_mesh.removed_points.append(cut_mesh.points[face[2]])
            cut_mesh.removed_points.append(cut_mesh.points[face[3]])

    cut_mesh.faces = np.array(tempFaces)
    cut_mesh.removed_faces = np.array(cut_mesh.removed_faces)
    cut_mesh.removed_points = np.array(cut_mesh.removed_points)
    cut_mesh._create_vtp_mesh()

    cut_mesh.files_base = os.path.splitext(cut_mesh.file_name)[0] + '_cut_mesh_bemio_output'

    print('Cut mesh in direction [' + str(plane_ind) + '] in direction [' + str(cut_dir) + '] at the location [' + str(plane_loc) + ']')

    return cut_mesh
