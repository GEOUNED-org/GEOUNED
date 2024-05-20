#
#  definition of GEOUNED objects to release memory
#
#  GEOUNED SurfacesGU, SolidsGU, PlaneGU, etc.. objects are created because FreeCAD create a new object
#  each time an attribute of FreeCAD object is called. This leads to code crash with memory failure
#  when attribues are call large amount of times. Like it is in this code.

import logging
import math

import Part

from .basic_functions_part1 import is_same_value
from .basic_functions_part2 import is_same_torus

logger = logging.getLogger("general_logger")


# SURFACES
class SurfacesGu(object):
    """GEOUNED surface class"""

    def __init__(self, face):
        self.face = face
        self.Surface = self.face.Surface
        self.type = str(self.Surface)

    def __str__(self):
        """str(Surface) is done for the classification of the surface.
        Surface_GU saves this information in self.type"""
        return self.type


class PlaneGu(SurfacesGu):
    """GEOUNED Plane Class"""

    def __init__(self, face, plane3Pts=False):
        SurfacesGu.__init__(self, face)
        self.Axis = face.Surface.Axis
        self.Position = face.Surface.Position
        self.pointDef = plane3Pts
        if plane3Pts:
            self.Points = tuple(v.Point for v in face.Vertexes)
            d1 = self.Points[0] - self.Points[1]
            d2 = self.Points[0] - self.Points[2]
            d3 = self.Points[1] - self.Points[2]
            self.dim1 = max(d1.Length, d2.Length, d3.Length)
            self.dim2 = min(d1.Length, d2.Length, d3.Length)
        else:
            self.dim1 = face.ParameterRange[1] - face.ParameterRange[0]
            self.dim2 = face.ParameterRange[3] - face.ParameterRange[2]


class CylinderGu(SurfacesGu):
    """GEOUNED Cylinder Class"""

    def __init__(self, face):
        SurfacesGu.__init__(self, face)
        self.Axis = face.Surface.Axis
        self.Radius = face.Surface.Radius
        self.Center = face.Surface.Center
        self.dimL = face.ParameterRange[3] - face.ParameterRange[2]


class ConeGu(SurfacesGu):
    """GEOUNED Cone Class"""

    def __init__(self, face):
        SurfacesGu.__init__(self, face)
        self.Axis = face.Surface.Axis
        self.Apex = face.Surface.Apex
        self.SemiAngle = face.Surface.SemiAngle
        self.dimL = face.ParameterRange[3] - face.ParameterRange[2]
        self.dimR = face.Surface.Radius
        self.Radius = face.Surface.Radius


class SphereGu(SurfacesGu):
    """GEOUNED Sphere Class"""

    def __init__(self, face):
        SurfacesGu.__init__(self, face)
        self.type = self.type[0:6]
        self.Center = face.Surface.Center
        self.Radius = face.Surface.Radius


class TorusGu(SurfacesGu):
    """GEOUNED Torus Class"""

    def __init__(self, face):
        SurfacesGu.__init__(self, face)
        self.Center = face.Surface.Center
        self.Axis = face.Surface.Axis
        self.MajorRadius = face.Surface.MajorRadius
        self.MinorRadius = face.Surface.MinorRadius


class SolidGu:
    """GEOUNED Solid Class"""

    def __init__(self, solid, tolerances, plane3Pts=False):
        self.solid = solid
        faces = define_list_face_gu(solid.Faces, plane3Pts)
        self.Faces = faces
        self.tolerances = tolerances
        self.Solids = solid.Solids
        self.BoundBox = solid.BoundBox
        self.Edges = solid.Edges
        self.TorusVParams = {}
        self.TorusUParams = {}

        toroidIndex = []
        for i, face in enumerate(self.Faces):
            if str(face.Surface) != "<Toroid object>":
                continue
            toroidIndex.append(i)

        if len(toroidIndex) != 0:
            tFaces = self.same_torus_surf(toroidIndex)
            for i, tSet in enumerate(tFaces):
                URange = self.merge_periodic_uv("U", tSet)
                VRange = self.merge_periodic_uv("V", tSet)
                for t in tSet:
                    self.TorusVParams[t] = (i, VRange)
                    self.TorusUParams[t] = (i, URange)

    def same_torus_surf(self, torusList):
        """group as a single face all the neighbour faces of the same torus"""
        sameTorusFace = []
        temp = torusList[:]
        while len(temp) > 0:
            i = temp[0]
            current = [i]
            for j in temp[1:]:
                if is_same_torus(
                    self.Faces[i].Surface,
                    self.Faces[j].Surface,
                    dtol=self.tolerances.tor_distance,
                    atol=self.tolerances.tor_angle,
                    rel_tol=self.tolerances.relativeTol,
                ):
                    current.append(j)
            for c in current:
                temp.remove(c)
            sameTorusFace.append(current)

        return self.separate_surfaces(sameTorusFace)

    def separate_surfaces(self, faceList):
        """group all faces in faceList forming a continuous surface"""
        sameSurfaces = []
        for tset in faceList:
            temp = tset[:]
            while len(temp) > 0:
                i = 0
                current = [temp[0]]
                removeList = [temp[0]]
                while len(temp) > 0 and i < len(current):
                    for tindex in temp:
                        if self.Faces[current[i]].distToShape(self.Faces[tindex])[0] < self.tolerances.distance:
                            if tindex not in current:
                                current.append(tindex)
                                removeList.append(tindex)
                    i += 1
                    for c in removeList:
                        temp.remove(c)
                    removeList = []

                sameSurfaces.append(current)
        return sameSurfaces

    # TODO check if this function is used as it appears to be nut used in the code
    def merge_no_periodic_uv(self, parameter, faceList):
        if parameter == "U":
            i1 = 0
            i2 = 2
        elif parameter == "V":
            i1 = 2
            i2 = 4

        v_min, v_max = self.Faces[faceList[0]].ParameterRange[i1:i2]
        for face in faceList[1:]:
            V0, V1 = self.Faces[face].ParameterRange[i1:i2]
            v_min = min(v_min, V0)
            v_max = max(v_max, V1)
        mergedParams = (False, (v_min, v_max))

        return mergedParams

    def merge_periodic_uv(self, parameter, faceList):
        two_pi = 2.0 * math.pi
        if parameter == "U":
            i1 = 0
            i2 = 2
        elif parameter == "V":
            i1 = 2
            i2 = 4

        params = []
        arcLength = 0.0
        for face in faceList:
            V0, V1 = self.Faces[face].ParameterRange[i1:i2]
            arcLength += V1 - V0
            params.append((V0, V1))

        params.sort()
        V0 = params[0][0]
        V1 = params[-1][1]
        if arcLength >= two_pi * (1.0 - self.tolerances.relativePrecision):
            mergedParams = (True, (V0, V0 + two_pi))
        else:
            if is_same_value(V0, 0.0, self.tolerances.relativePrecision) and is_same_value(
                V1, two_pi, self.tolerances.relativePrecision
            ):
                for i in range(len(params) - 1):
                    if not is_same_value(
                        params[i][1],
                        params[i + 1][0],
                        self.tolerances.relativePrecision,
                    ):
                        break
                v_min = params[i + 1][0] - two_pi
                v_max = params[i][1]
            else:
                v_min = params[0][0]
                v_max = params[-1][1]
            mergedParams = (False, (v_min, v_max))

        return mergedParams


# FACES
class FaceGu(object):
    """GEOUNED Face Class"""

    def __init__(self, face, Plane3Pts=False):
        # GEOUNED based atributes
        self.__face__ = face
        self.Surface = define_surface(face, Plane3Pts)  # Define the appropiate GU Surface of the face

        # FreeCAD based Atributes
        self.Area = face.Area
        self.CenterOfMass = face.CenterOfMass
        self.ParameterRange = face.ParameterRange
        self.Orientation = face.Orientation
        self.OuterWire = face.OuterWire
        self.Edges = face.Edges
        self.Vertexes = face.Vertexes
        return

    def tessellate(self, val, reset=False):
        res = self.__face__.tessellate(val, reset)
        return res

    def getUVNodes(self):
        return self.__face__.getUVNodes()

    def isEqual(self, face):
        return self.__face__.isEqual(face.__face__)

    def valueAt(self, u, v):
        return self.__face__.valueAt(u, v)

    def distToShape(self, shape):
        shape1 = self.__face__
        if type(shape) is Part.Shape:
            shape2 = shape
        else:
            shape2 = shape.__face__

        if shape1 is shape2:
            return (0,)
        else:
            try:
                dist2Shape = shape1.distToShape(shape2)
            except:
                dist2Shape = shape2.distToShape(shape1)
            return dist2Shape


# Aux functions
def define_list_face_gu(face_list, plane3Pts=False):
    """Return the list of the  corresponding Face_GU  object of a FaceList"""
    return tuple(FaceGu(face, plane3Pts) for face in face_list)


def define_surface(face, plane3Pts):
    kind_surf = str(face.Surface)
    if kind_surf == "<Plane object>":
        Surf_GU = PlaneGu(face, plane3Pts)
    elif kind_surf == "<Cylinder object>":
        Surf_GU = CylinderGu(face)
    elif kind_surf == "<Cone object>":
        Surf_GU = ConeGu(face)
    elif kind_surf[0:6] == "Sphere":
        Surf_GU = SphereGu(face)
    elif kind_surf == "<Toroid object>":
        Surf_GU = TorusGu(face)
    else:
        logger.info(f"bad Surface type {kind_surf}")
        Surf_GU = None
    return Surf_GU
