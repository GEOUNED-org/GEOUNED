#
#  definition of GEOUNED objects to release memory
#
#  GEOUNED SurfacesGU, SolidsGU, PlaneGU, etc.. objects are created because FreeCAD create a new object
#  each time an attribute of FreeCAD object is called. This leads to code crash with memory failure
#  when attribues are call large amount of times. Like it is in this code.

import math
import Part

from ..Utils.Options.Classes import Tolerances as tol
from .BasicFunctions_part1 import isSameValue
from .BasicFunctions_part2 import isSameTorus


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

    def __init__(self, solid, plane3Pts=False):
        self.solid = solid
        faces = DefineListFace_GU(solid.Faces, plane3Pts)
        self.Faces = faces
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
            tFaces = self.__sameTorusSurf__(toroidIndex)
            for i, tSet in enumerate(tFaces):
                URange = self.__mergePeriodicUV__("U", tSet)
                VRange = self.__mergePeriodicUV__("V", tSet)
                for t in tSet:
                    self.TorusVParams[t] = (i, VRange)
                    self.TorusUParams[t] = (i, URange)

    def __sameTorusSurf__(self, torusList):
        """group as a single face all the neighbour faces of the same torus"""
        sameTorusFace = []
        temp = torusList[:]
        while len(temp) > 0:
            i = temp[0]
            current = [i]
            for j in temp[1:]:
                if isSameTorus(
                    self.Faces[i].Surface,
                    self.Faces[j].Surface,
                    dtol=tol.tor_distance,
                    atol=tol.tor_angle,
                    relTol=tol.relativeTol,
                ):
                    current.append(j)
            for c in current:
                temp.remove(c)
            sameTorusFace.append(current)

        return self.__separateSurfaces__(sameTorusFace)

    def __separateSurfaces__(self, faceList):
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
                        if (
                            self.Faces[current[i]].distToShape(self.Faces[tindex])[0]
                            < tol.distance
                        ):
                            if tindex not in current:
                                current.append(tindex)
                                removeList.append(tindex)
                    i += 1
                    for c in removeList:
                        temp.remove(c)
                    removeList = []

                sameSurfaces.append(current)
        return sameSurfaces

    def __mergeNoPeriodicUV__(self, parameter, faceList):
        if parameter == "U":
            i1 = 0
            i2 = 2
        elif parameter == "V":
            i1 = 2
            i2 = 4

        Vmin, Vmax = self.Faces[faceList[0]].ParameterRange[i1:i2]
        for face in faceList[1:]:
            V0, V1 = self.Faces[face].ParameterRange[i1:i2]
            Vmin = min(Vmin, V0)
            Vmax = max(Vmax, V1)
        mergedParams = (False, (Vmin, Vmax))

        return mergedParams

    def __mergePeriodicUV__(self, parameter, faceList):
        twoPi = 2.0 * math.pi
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
        if arcLength >= twoPi * (1.0 - tol.relativePrecision):
            mergedParams = (True, (V0, V0 + twoPi))
        else:
            if isSameValue(V0, 0.0, tol.relativePrecision) and isSameValue(
                V1, twoPi, tol.relativePrecision
            ):
                for i in range(len(params) - 1):
                    if not isSameValue(
                        params[i][1], params[i + 1][0], tol.relativePrecision
                    ):
                        break
                Vmin = params[i + 1][0] - twoPi
                Vmax = params[i][1]
            else:
                Vmin = params[0][0]
                Vmax = params[-1][1]
            mergedParams = (False, (Vmin, Vmax))

        return mergedParams


# FACES
class FaceGu(object):
    """GEOUNED Face Class"""

    def __init__(self, face, Plane3Pts=False):
        # GEOUNED based atributes
        self.__face__ = face
        self.Surface = DefineSurface(
            face, Plane3Pts
        )  # Define the appropiate GU Surface of the face

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
def DefineListFace_GU(face_list, plane3Pts=False):
    """Return the list of the  corresponding Face_GU  object of a FaceList"""
    return tuple(FaceGu(face, plane3Pts) for face in face_list)


def DefineSurface(face, plane3Pts):
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
        print("bad Surface type", kind_surf)
        Surf_GU = None
    return Surf_GU


def ListSurfaces(Surfaces):
    Faces = []
    for elem in Surfaces:
        Faces.extend(DefineSurface(face) for face in elem)
    return Faces
