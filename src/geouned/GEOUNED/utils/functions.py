#
# Set of useful functions used in different parts of the code
#
import logging
import math

import BOPTools.SplitAPI
import FreeCAD
import numpy as np
import Part

logger = logging.getLogger("general_logger")

from .basic_functions_part1 import (
    ConeParams,
    CylinderParams,
    Plane3PtsParams,
    PlaneParams,
    SphereParams,
    TorusParams,
    is_parallel,
)
from . import basic_functions_part2 as BF


def get_box(comp, options):
    bb = FreeCAD.BoundBox(comp.BoundBox)
    bb.enlarge(options.enlargeBox)
    xMin, yMin, zMin = bb.XMin, bb.YMin, bb.ZMin
    xLength, yLength, zLength = bb.XLength, bb.YLength, bb.ZLength

    return Part.makeBox(
        xLength,
        yLength,
        zLength,
        FreeCAD.Vector(xMin, yMin, zMin),
        FreeCAD.Vector(0, 0, 1),
    )


class GeounedSolid:

    def __init__(self, id, comsolid=None):
        refine = True
        if not comsolid:
            self.Solids = None
            self.Volume = None
            self.BoundBox = None
        elif type(comsolid) is list:
            self.Solids = comsolid
            vol = 0
            self.BoundBox = FreeCAD.BoundBox()
            for s in comsolid:
                vol += s.Volume
                self.BoundBox.add(s.BoundBox)
            self.Volume = vol
        else:
            if refine:
                try:
                    self.Solids = comsolid.removeSplitter().Solids
                except:
                    self.Solids = comsolid.Solids

                for s in self.Solids:
                    if s.Volume < 0:
                        s.reverse()
            else:
                self.Solids = comsolid.Solids
            self.Volume = comsolid.Volume
            self.BoundBox = comsolid.BoundBox

        self.__id__ = id
        self.label = None
        self.Definition = []
        self.Faces = []
        self.Comments = ""
        self.Density = 0
        self.Dilution = 1
        self.Material = 0
        self.Surfaces = []
        self.Rho = None
        self.MatInfo = None
        self.CellType = "solid"  # other types : 'void', 'enclosure', 'fill'
        self.Universe = 0
        self.Void = False
        self.IsEnclosure = False
        self.EnclosureID = None
        self.ParentEnclosureID = None
        self.SonEnclosures = []
        self.enclosure_list = None
        self.CADSolid = None
        self.UniverseBox = None
        self.NullCell = True

    def update_solids(self, solidList):
        self.Solids = solidList
        vol = 0
        self.BoundBox = FreeCAD.BoundBox()
        for s in solidList:
            vol += s.Volume
            self.BoundBox.add(s.BoundBox)

    def set_cad_solid(self):
        self.CADSolid = Part.makeCompound(self.Solids)
        self.Volume = self.CADSolid.Volume
        self.BoundBox = self.CADSolid.BoundBox

    def set_definition(self, definition, simplify=False):

        if definition is None:
            self.NullCell = True
            return

        self.NullCell = False
        self.Definition = definition

        if not self.Void:
            if not self.Definition.elements:
                self.NullCell = True
                return

        self.Surfaces = tuple(self.Definition.get_surfaces_numbers())

    def set_faces(self, faces):
        self.Faces = faces

    def set_comments(self, comments):
        self.Comments = comments

    def set_material(self, material, rho=None, info=None):

        self.Material = material
        self.Rho = rho
        self.MatInfo = info

        if rho is not None:
            self.Density = self.Dilution * rho
        else:
            if material != 0:
                self.Density = None

    def set_dilution(self, dilution):
        self.Dilution = dilution
        if self.Rho is not None:
            self.Density = self.Rho * dilution

    def check_intersection(self, solid, dtolerance=1.0e-6, vtolerance=1e-10):
        """Check if solid intersect with current solid.
        return : -2 solid fully embedded in self.CADSolid ;
                 -1 self.CADSolid fully embedded in solid ;
                  0 self.CADSolid intersect solid ;
                  1 self.CADSolid and solid fully disjoint"""

        dist = 1e12
        for sol1 in self.CADSolid.Solids:
            for sol2 in solid.Solids:
                try:
                    distShape = sol1.distToShape(sol2)[0]
                except:
                    logger.info("Failed solid1.distToshape(solid2), try with inverted solids")
                    distShape = sol2.distToShape(sol1)[0]
                    logger.info(f"inverted disToShape OK {distShape}")
                dist = min(dist, distShape)
                if dist == 0:
                    break

        if dist > dtolerance:
            return 1

        common = self.CADSolid.common(solid)
        if abs(common.Volume) < vtolerance:
            return 1
        if abs(self.CADSolid.Volume - common.Volume) / common.Volume < vtolerance:
            return -1
        elif abs(solid.Volume - common.Volume) / common.Volume < vtolerance:
            return -2
        else:
            return 0


class GeounedSurface:

    def __init__(self, params, boundBox, Face=None):

        self.Index = 0
        self.__boundBox__ = boundBox
        if params[0] == "Plane":
            self.Type = "Plane"
            self.Surf = PlaneParams(params[1])  # plane point defined as the shortest distance to origin
        elif params[0] == "Plane3Pts":
            self.Type = "Plane"
            self.Surf = Plane3PtsParams(params[1])  # plane point defined with 3 points
        elif params[0] == "Cylinder":
            self.Type = params[0]
            self.Surf = CylinderParams(params[1])
        elif params[0] == "Cone":
            self.Type = params[0]
            self.Surf = ConeParams(params[1])
        elif params[0] == "Sphere":
            self.Type = params[0]
            self.Surf = SphereParams(params[1])
        elif params[0] == "Torus":
            self.Type = params[0]
            self.Surf = TorusParams(params[1])

        self.shape = Face

        if type(Face) is str:
            if Face == "Build":
                self.build_surface()

        if self.shape == "Build":
            raise ValueError(f"stop {params} {boundBox}")
        return

    def build_surface(self):

        if self.Type == "Plane":

            normal = self.Surf.Axis
            p0 = normal.dot(self.Surf.Position)
            Box = FreeCAD.BoundBox(self.__boundBox__)
            Box.enlarge(10)

            pointEdge = []
            for i in range(12):
                edge = Box.getEdge(i)
                p1 = normal.dot(edge[0])
                p2 = normal.dot(edge[1])
                d0 = p0 - p1
                d1 = p2 - p1
                if d1 != 0:
                    a = d0 / d1
                    if a >= 0 and a <= 1:
                        pointEdge.append(edge[0] + a * (edge[1] - edge[0]))

            if len(pointEdge) == 0:
                self.shape = None  # Plane does not cross box
                return

            s = FreeCAD.Vector((0, 0, 0))
            for v in pointEdge:
                s = s + v
            s = s / len(pointEdge)

            vtxvec = []
            for v in pointEdge:
                vtxvec.append(v - s)

            X0 = vtxvec[0]
            Y0 = normal.cross(X0)

            orden = []
            for i, v in enumerate(vtxvec):
                phi = np.arctan2(v.dot(Y0), v.dot(X0))
                orden.append((phi, i))
            orden.sort()

            self.shape = Part.Face(Part.makePolygon([pointEdge[p[1]] for p in orden], True))

            return

        elif self.Type == "Cylinder":
            dir = self.Surf.Axis
            pnt = self.Surf.Center
            rad = self.Surf.Radius

            dmin = dir.dot((self.__boundBox__.getPoint(0) - pnt))
            dmax = dmin
            for i in range(1, 8):
                d = dir.dot((self.__boundBox__.getPoint(i) - pnt))
                dmin = min(d, dmin)
                dmax = max(d, dmax)

            pnt = pnt + (dmin - 5) * dir
            length = dmax - dmin + 10
            self.shape = Part.makeCylinder(rad, length, pnt, dir, 360.0)
            return

        elif self.Type == "Cone":
            #         dir  = self.Surf.Axis
            #         apex = self.Surf.Apex
            #         tan = math.tan(self.Surf.SemiAngle)
            #         rad = tan*diag # New radius far away to cover all geometry

            axis = self.Surf.Axis
            apex = self.Surf.Apex
            tan = math.tan(self.Surf.SemiAngle)

            dmax = axis.dot((self.__boundBox__.getPoint(0) - apex))
            for i in range(1, 8):
                d = axis.dot((self.__boundBox__.getPoint(i) - apex))
                dmax = max(d, dmax)

            if dmax > 0:
                length = dmax + 10
                rad = tan * length
                self.shape = Part.makeCone(0.0, rad, length, apex, axis, 360.0)
            else:
                self.shape = None
            return

        elif self.Type == "Sphere":
            rad = self.Surf.Radius
            pnt = self.Surf.Center
            self.shape = Part.makeSphere(rad, pnt).Faces[0]
            return

        elif self.Type == "Torus":
            axis = self.Surf.Axis
            center = self.Surf.Center
            majorR = self.Surf.MajorRadius
            minorR = self.Surf.MinorRadius

            torus = Part.makeTorus(majorR, minorR, center, axis)
            self.shape = torus.Faces[0]
            return
        else:
            logger.error(f"Type {self.Type} is not defined")
            return


class SurfacesDict(dict):
    def __init__(self, surfaces=None, offset=0):
        self.IndexOffset = offset
        surfname = ["PX", "PY", "PZ", "P", "Cyl", "Cone", "Sph", "Tor"]
        for name in surfname:
            self[name] = []

        self.__surfIndex__ = dict()

        if surfaces is not None:
            for key in surfaces.keys():
                self[key] = surfaces[key][:]
                self.__surfIndex__[key] = surfaces.__surfIndex__[key][:]
            self.surfaceNumber = surfaces.surfaceNumber
            self.__last_obj__ = (surfaces.__last_obj__[0], surfaces.__last_obj__[1])
        else:
            self.surfaceNumber = 0
            self.__last_obj__ = ("", -1)
            for key in surfname:
                self.__surfIndex__[key] = []
        return

    def __str__(self):
        for key in self.keys():
            logger.info(f"{key}, {self[key]}")
        return ""

    def get_surface(self, index):

        lastKey = self.__last_obj__[0]
        lastInd = self.__last_obj__[1]
        if lastKey != "":
            if len(self[lastKey]) > 0:
                if self[lastKey][lastInd].Index == index:
                    return self[lastKey][lastInd]

        for key, values in self.__surfIndex__.items():
            if index not in values:
                continue
            i = values.index(index)
            self.__last_obj__ = (key, i)
            return self[key][i]

        logger.info(f"Index {index} not found in Surfaces")
        return None

    def del_surface(self, index):
        self.get_surface(index)
        self.__surfIndex__[self.__last_obj__[0]].remove(index)
        del self[self.__last_obj__[0]][self.__last_obj__[1]]
        return

    def extend(self, surface, options, tolerances, numeric_format):
        for Pkey in ["PX", "PY", "PZ", "P"]:
            for s in surface[Pkey]:
                self.add_plane(s, options, tolerances, numeric_format, False)
        for s in surface["Cyl"]:
            self.add_cylinder(s, options, tolerances, numeric_format, False)
        for s in surface["Cone"]:
            self.add_cone(s, tolerances)
        for s in surface["Sph"]:
            self.add_sphere(s, tolerances)
        for s in surface["Tor"]:
            self.add_torus(s, tolerances)

    def add_plane(self, plane, options, tolerances, numeric_format, fuzzy):
        ex = FreeCAD.Vector(1, 0, 0)
        ey = FreeCAD.Vector(0, 1, 0)
        ez = FreeCAD.Vector(0, 0, 1)

        if is_parallel(plane.Surf.Axis, ex, tolerances.pln_angle):
            add_plane = True
            for i, p in enumerate(self["PX"]):
                if BF.is_same_plane(
                    plane.Surf,
                    p.Surf,
                    options=options,
                    tolerances=tolerances,
                    numeric_format=numeric_format,
                    fuzzy=(fuzzy, p.Index),
                ):
                    add_plane = False
                    index = p.Index
                    self.__last_obj__ = ("PX", i)
                    break
            if add_plane:
                self.surfaceNumber += 1
                plane.Index = self.surfaceNumber + self.IndexOffset
                self.__last_obj__ = ("PX", len(self["PX"]))
                self["PX"].append(plane)
                self.__surfIndex__["PX"].append(plane.Index)

        elif is_parallel(plane.Surf.Axis, ey, tolerances.pln_angle):
            add_plane = True
            for i, p in enumerate(self["PY"]):
                if BF.is_same_plane(
                    plane.Surf,
                    p.Surf,
                    options=options,
                    tolerances=tolerances,
                    numeric_format=numeric_format,
                    fuzzy=(fuzzy, p.Index),
                ):
                    add_plane = False
                    index = p.Index
                    self.__last_obj__ = ("PY", i)
                    break
            if add_plane:
                self.surfaceNumber += 1
                plane.Index = self.surfaceNumber + self.IndexOffset
                self.__last_obj__ = ("PY", len(self["PY"]))
                self["PY"].append(plane)
                self.__surfIndex__["PY"].append(plane.Index)

        elif is_parallel(plane.Surf.Axis, ez, tolerances.pln_angle):
            add_plane = True
            for i, p in enumerate(self["PZ"]):
                if BF.is_same_plane(
                    plane.Surf,
                    p.Surf,
                    options=options,
                    tolerances=tolerances,
                    numeric_format=numeric_format,
                    fuzzy=(fuzzy, p.Index),
                ):
                    add_plane = False
                    index = p.Index
                    self.__last_obj__ = ("PZ", i)
                    break
            if add_plane:
                self.surfaceNumber += 1
                plane.Index = self.surfaceNumber + self.IndexOffset
                self.__last_obj__ = ("PZ", len(self["PZ"]))
                self["PZ"].append(plane)
                self.__surfIndex__["PZ"].append(plane.Index)

        else:
            add_plane = True
            for i, p in enumerate(self["P"]):
                if BF.is_same_plane(
                    plane.Surf,
                    p.Surf,
                    options=options,
                    tolerances=tolerances,
                    numeric_format=numeric_format,
                    fuzzy=(fuzzy, p.Index),
                ):
                    add_plane = False
                    index = p.Index
                    self.__last_obj__ = ("P", i)
                    break
            if add_plane:
                self.surfaceNumber += 1
                plane.Index = self.surfaceNumber + self.IndexOffset
                self.__last_obj__ = ("P", len(self["P"]))
                self["P"].append(plane)
                self.__surfIndex__["P"].append(plane.Index)

        if add_plane:
            return plane.Index, False
        else:
            return index, True

    def add_cylinder(self, cyl, options, tolerances, numeric_format, fuzzy=False):
        addCyl = True
        for i, c in enumerate(self["Cyl"]):
            if BF.is_same_cylinder(
                cyl.Surf,
                c.Surf,
                options=options,
                tolerances=tolerances,
                numeric_format=numeric_format,
                fuzzy=(fuzzy, c.Index),
            ):
                addCyl = False
                index = c.Index
                self.__last_obj__ = ("Cyl", i)
                break

        if addCyl:
            self.surfaceNumber += 1
            cyl.Index = self.surfaceNumber + self.IndexOffset
            self.__last_obj__ = ("Cyl", len(self["Cyl"]))
            self["Cyl"].append(cyl)
            self.__surfIndex__["Cyl"].append(cyl.Index)
            return cyl.Index, False
        else:
            return index, True

    def add_cone(self, cone, tolerances):
        cone_added = True
        for i, c in enumerate(self["Cone"]):
            if BF.is_same_cone(
                cone.Surf,
                c.Surf,
                dtol=tolerances.kne_distance,
                atol=tolerances.kne_angle,
                rel_tol=tolerances.relativeTol,
            ):
                cone_added = False
                index = c.Index
                self.__last_obj__ = ("Cone", i)
                break
        if cone_added:
            self.surfaceNumber += 1
            cone.Index = self.surfaceNumber + self.IndexOffset
            self.__last_obj__ = ("Cone", len(self["Cone"]))
            self["Cone"].append(cone)
            self.__surfIndex__["Cone"].append(cone.Index)
            return cone.Index, False
        else:
            return index, True

    def add_sphere(self, sph, tolerances):
        sphere_added = True
        for i, s in enumerate(self["Sph"]):
            if BF.is_same_sphere(
                sph.Surf,
                s.Surf,
                tolerances.sph_distance,
                rel_tol=tolerances.relativeTol,
            ):
                sphere_added = False
                index = s.Index
                self.__last_obj__ = ("Sph", i)
                break
        if sphere_added:
            self.surfaceNumber += 1
            sph.Index = self.surfaceNumber + self.IndexOffset
            self.__last_obj__ = ("Sph", len(self["Sph"]))
            self["Sph"].append(sph)
            self.__surfIndex__["Sph"].append(sph.Index)
            return sph.Index, False
        else:
            return index, True

    def add_torus(self, tor, tolerances):
        add_torus = True
        for i, s in enumerate(self["Tor"]):
            if BF.is_same_torus(
                tor.Surf,
                s.Surf,
                dtol=tolerances.tor_distance,
                atol=tolerances.tor_angle,
                rel_tol=tolerances.relativeTol,
            ):
                add_torus = False
                index = s.Index
                self.__last_obj__ = ("Tor", i)
                break
        if add_torus:
            self.surfaceNumber += 1
            tor.Index = self.surfaceNumber + self.IndexOffset
            self.__last_obj__ = ("Tor", len(self["Tor"]))
            self["Tor"].append(tor)
            self.__surfIndex__["Tor"].append(tor.Index)
            return tor.Index, False
        else:
            return index, True


def split_bop(solid, tools, tolerance, options, scale=0.1):

    if tolerance >= 0.1:
        compSolid = BOPTools.SplitAPI.slice(solid, tools, "Split", tolerance=tolerance)

    elif tolerance < 1e-12:
        if options.scaleUp:
            tol = 1e-13 if options.splitTolerance == 0 else options.splitTolerance
            compSolid = split_bop(solid, tools, tol / scale, options, 1.0 / scale)
        else:
            compSolid = BOPTools.SplitAPI.slice(solid, tools, "Split", tolerance=tolerance)

    else:
        try:
            compSolid = BOPTools.SplitAPI.slice(solid, tools, "Split", tolerance=tolerance)
        except:
            compSolid = split_bop(solid, tools, tolerance * scale, options, scale)

    return compSolid
