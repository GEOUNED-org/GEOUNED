#   Conversion to MCNP v0.0
#   Only one solid and planar surfaces
#

import logging
import math
from collections import OrderedDict

import FreeCAD
import Part

from ..conversion import cell_definition as CD
from ..utils import functions as UF
from ..utils import geometry_gu as GU
from ..utils.basic_functions_part1 import (
    is_in_line,
    is_in_plane,
    is_parallel,
    is_same_value,
)
from ..utils.basic_functions_part2 import is_duplicate_in_list

logger = logging.getLogger("general_logger")

twoPi = math.pi * 2


def split_full_cylinder(solid, options, tolerances, numeric_format):
    explode = []
    bases = [solid]
    while True:
        new_bases = []
        for base in bases:
            cut_solids = cut_full_cylinder(base, options, tolerances, numeric_format)
            if len(cut_solids) == 1:
                explode.extend(cut_solids)
            else:
                new_bases.extend(cut_solids)
        if len(new_bases) == 0:
            break
        else:
            bases = new_bases

    return Part.makeCompound(explode)


def cut_full_cylinder(solid, options, tolerances, numeric_format):
    solid_gu = GU.SolidGu(solid, tolerances=tolerances)
    surfaces = UF.SurfacesDict()
    flag_inv = CD.is_inverted(solid_gu.solid)
    universe_box = solid.BoundBox

    for face in solid_gu.Faces:
        surf = str(face.Surface)
        if surf == "<Cylinder object>":
            if flag_inv:
                orient = "Reversed" if face.Orientation == "Forward" else "Forward"
            else:
                orient = face.Orientation

            u1, u2, v1, v2 = face.ParameterRange
            angle = abs(u2 - u1)

            # closed convex cylinder
            if abs(angle % twoPi) < 1e-2 and orient == "Forward":
                dir = face.Surface.Axis
                orig = face.Surface.Center
                rad = face.Surface.Radius
                dim_l = face.ParameterRange[3] - face.ParameterRange[2]
                cylinder = UF.GeounedSurface(("Cylinder", (orig, dir, rad, dim_l)), universe_box)
                cylinder.build_surface()
                surfaces.add_cylinder(cylinder, options, tolerances, numeric_format, False)

                # add planes if cylinder axis is cut by a plane (plane quasi perpedicular to axis)
                for p in cyl_bound_planes(face, universe_box):
                    p.build_surface()
                    surfaces.add_plane(p, options, tolerances, numeric_format, False)
                break

    planes = []
    for P in ("PX", "PY", "PZ", "P"):
        planes.extend(surfaces[P])
    planes = sort_planes(planes, True)

    if len(planes) == 0:
        return [solid]
    if len(planes[-1]) < options.nPlaneReverse:
        planes.reverse()

    cut = False
    for pp in planes:
        # cut with more external parallel planes
        pp[0].build_surface()
        if len(pp) != 1:
            pp[-1].build_surface()
            tools = (pp[0].shape, pp[-1].shape)
        else:
            tools = (pp[0].shape,)

        try:
            comsolid = UF.split_bop(solid, tools, options.splitTolerance, options)
        except:
            comsolid = solid
        if len(comsolid.Solids) > 1:
            out_solid = comsolid.Solids
            cut = True
            break

    if cut:
        return out_solid

    tool = (surfaces["Cyl"][0].shape,)
    try:
        comsolid = UF.split_bop(solid, tool, options.splitTolerance, options)
    except:
        comsolid = solid

    if len(comsolid.Solids) > 1:
        out_solid = comsolid.Solids
    else:
        out_solid = [solid]

    return out_solid


def gen_plane(pos, normal, diag):
    plane = Part.makePlane(diag, diag, pos, normal)
    vec_on_plane = plane.Vertexes[3].Point.sub(plane.Vertexes[0].Point)
    new_pos = plane.Vertexes[0].Point.sub(vec_on_plane)
    plane_center = Part.makePlane(2.0 * diag, 2.0 * diag, new_pos, normal)
    return plane_center


def cyl_bound_planes(face, boundBox):
    Edges = face.OuterWire.Edges
    planes = []
    for e in Edges:
        try:
            curve = str(e.Curve)
        except:
            curve = "none"

        if curve[0:6] == "Circle":
            dir = e.Curve.Axis
            center = e.Curve.Center
            dim1 = e.Curve.Radius
            dim2 = e.Curve.Radius
            plane = UF.GeounedSurface(("Plane", (center, dir, dim1, dim2)), boundBox)
            planes.append(plane)

        elif curve == "<Ellipse object>":
            dir = e.Curve.Axis
            center = e.Curve.Center
            dim1 = e.Curve.MinorRadius
            dim2 = e.Curve.MajorRadius
            plane = UF.GeounedSurface(("Plane", (center, dir, dim1, dim2)), boundBox)
            planes.append(plane)

    return planes


def torus_bound_planes(face, boundBox, tolerances):
    params = face.ParameterRange
    planes = []
    if is_same_value(params[1] - params[0], twoPi, tolerances.value):
        return planes

    Edges = face.OuterWire.Edges

    for e in Edges:
        try:
            curve = str(e.Curve)
        except:
            curve = "none"

        if curve[0:6] == "Circle":
            dir = e.Curve.Axis
            if not is_parallel(dir, face.Surface.Axis, tolerances.angle):
                center = e.Curve.Center
                dim1 = e.Curve.Radius
                dim2 = e.Curve.Radius
                plane = UF.GeounedSurface(("Plane", (center, dir, dim1, dim2)), boundBox)
                planes.append(plane)

        elif curve == "<Ellipse object>":
            dir = e.Curve.Axis
            center = e.Curve.Center
            dim1 = e.Curve.MinorRadius
            dim2 = e.Curve.MajorRadius
            plane = UF.GeounedSurface(("Plane", (center, dir, dim1, dim2)), boundBox)
            planes.append(plane)

        elif curve == "<BSplineCurve object>":
            planeParams = plane_spline_curve(e, tolerances)
            if planeParams is not None:
                plane = UF.GeounedSurface(("Plane", planeParams), boundBox)
                planes.append(plane)

    return planes


def plane_spline_curve(edge, tolerances):

    normal = edge.derivative1At(0).cross(edge.derivative1At(0.5))
    normal.normalize()
    curve_2d = True
    for p in (0.25, 0.75, 1):
        # check if derivative orthogonal to curve normal vector
        if abs(normal.dot(edge.derivative1At(p))) > tolerances.value:
            curve_2d = False
            break

    r = edge.valueAt(0.25) - edge.valueAt(0.75)
    if curve_2d:
        return (edge.valueAt(0), normal, r.Length, r.Length)
    else:
        return None


def extract_surfaces(solid, kind, universe_box, options, tolerances, numeric_format, MakeObj):
    if kind == "All":
        fuzzy = True
        solid_parts = []
        for s in solid.Solids:
            solid_parts.append(GU.SolidGu(s, tolerances=tolerances))
    else:
        fuzzy = False
        if kind == "Plane3Pts":
            P3P = True
        else:
            P3P = False
        solid_parts = [GU.SolidGu(solid, tolerances=tolerances, plane3Pts=P3P)]

    surfaces = UF.SurfacesDict()

    for solid_GU in solid_parts:

        # Get the parameter of all faces of the solid
        # Add auxillary planes for Cylinders and Cones
        for face in solid_GU.Faces:
            surf = str(face.Surface)
            if surf == "<Plane object>" and kind in ["Planes", "All"]:
                normal = face.Surface.Axis
                pos = face.CenterOfMass
                dim1 = face.ParameterRange[1] - face.ParameterRange[0]
                dim2 = face.ParameterRange[3] - face.ParameterRange[2]
                plane = UF.GeounedSurface(("Plane", (pos, normal, dim1, dim2)), universe_box)
                if MakeObj:
                    plane.build_surface()
                surfaces.add_plane(plane, options, tolerances, numeric_format, fuzzy)

            elif surf == "<Cylinder object>":
                dir = face.Surface.Axis
                orig = face.Surface.Center
                rad = face.Surface.Radius
                dim_l = face.ParameterRange[3] - face.ParameterRange[2]
                if kind in ["Cyl", "All"]:
                    cylinder = UF.GeounedSurface(("Cylinder", (orig, dir, rad, dim_l)), universe_box)
                    if MakeObj:
                        cylinder.build_surface()
                    surfaces.add_cylinder(cylinder, options, tolerances, numeric_format, fuzzy)

                if kind in ["Planes", "All"]:
                    # add planes if cylinder axis is cut by a plane (plane quasi perpedicular to axis)
                    for p in cyl_bound_planes(face, universe_box):
                        if MakeObj:
                            p.build_surface()
                        surfaces.add_plane(p, options, tolerances, numeric_format, False)

            elif surf == "<Cone object>":
                dir = face.Surface.Axis
                apex = face.Surface.Apex
                half_angle = face.Surface.SemiAngle
                dim_l = face.ParameterRange[3] - face.ParameterRange[2]
                dimR = face.Surface.Radius
                if kind in ["Cone", "All"]:
                    cone = UF.GeounedSurface(("Cone", (apex, dir, half_angle, dim_l, dimR)), universe_box)
                    if MakeObj:
                        cone.build_surface()
                    surfaces.add_cone(cone, tolerances)

                if kind in ["Planes", "All"]:
                    for p in cyl_bound_planes(face, universe_box):
                        if MakeObj:
                            p.build_surface()
                        surfaces.add_plane(p, options, tolerances, numeric_format, False)

            elif surf[0:6] == "Sphere" and kind in ["Sph", "All"]:
                rad = face.Surface.Radius
                pnt = face.Surface.Center
                sphere = UF.GeounedSurface(("Sphere", (pnt, rad)), universe_box)
                if MakeObj:
                    sphere.build_surface()
                surfaces.add_sphere(sphere, tolerances)

                if kind in ["Planes", "All"]:
                    for p in cyl_bound_planes(face, universe_box):
                        if MakeObj:
                            p.build_surface()
                        surfaces.add_plane(p, options, tolerances, numeric_format, False)

            elif surf == "<Toroid object>":
                if kind in ["Tor", "All"]:
                    radMaj = face.Surface.MajorRadius
                    radMin = face.Surface.MinorRadius
                    center = face.Surface.Center
                    dir = face.Surface.Axis
                    torus = UF.GeounedSurface(("Torus", (center, dir, radMaj, radMin)), universe_box)
                    if MakeObj:
                        torus.build_surface()
                    surfaces.add_torus(torus, tolerances)

                if kind in ["Planes", "All"]:
                    for p in torus_bound_planes(face, universe_box, tolerances):
                        if MakeObj:
                            p.build_surface()
                        surfaces.add_plane(p, options, tolerances, numeric_format, False)

            elif surf == "<Plane object>" and kind == "Plane3Pts":
                pos = face.CenterOfMass
                normal = face.Surface.Axis
                dim1 = face.ParameterRange[1] - face.ParameterRange[0]
                dim2 = face.ParameterRange[3] - face.ParameterRange[2]
                points = tuple(v.Point for v in face.Vertexes)

                plane = UF.GeounedSurface(("Plane3Pts", (pos, normal, dim1, dim2, points)), universe_box)
                if MakeObj:
                    plane.build_surface()
                surfaces.add_plane(plane, options, tolerances, numeric_format, fuzzy)

    return surfaces


#   Check if to faces are joint
def contiguous_face(face1, face2, tolerances):
    return face1.distToShape(face2)[0] < tolerances.distance


def same_faces(Faces, tolerances):
    Connection = OrderedDict()
    if len(Faces) == 1:
        return []

    for i, face1 in enumerate(Faces):
        Couples = []
        if not Faces[i + 1 :]:
            continue
        for j, face2 in enumerate(Faces[i + 1 :]):
            if contiguous_face(face1, face2, tolerances):

                Couples.append(i + 1 + j)

        Connection[i] = Couples

    lista = Connection[0]
    Connection.popitem(0)

    if len(Connection) == 0:  # solo los elementos de la lista de la superficie 0
        return lista

    if not lista:  # ninguna face está conecta conecta con la superficie 0
        return lista

    for elem in Connection:
        if elem in lista:  # que la key esta en lista implica que sus dependencias estan
            lista.extend(Connection[elem])
        else:
            for elem2 in Connection[elem]:
                if elem2 in lista:  # si una de sus dependencias esta en lista lo esta la clave
                    lista.append(elem)

    return list(set(lista))


# Tolerance in this function are not the general once
# function should be reviewed
def gen_plane_cylinder(face, solid, tolerances):
    Surf = face.Surface
    rad = Surf.Radius
    if face.Area < 1e-2:
        return None

    UVNodes = []
    face_index_0 = [solid.Faces.index(face)]

    try:
        face.tessellate(0.1)
        UVNodes.append(face.getUVNodes())
    except RuntimeError:
        PR = face.ParameterRange
        UVNode1 = (PR[0], PR[2])
        UVNode2 = (PR[1], PR[3])
        UVNodes.append([UVNode1, UVNode2])

    for i, face2 in enumerate(solid.Faces):

        if str(face2.Surface) == "<Cylinder object>" and not (face2.isEqual(face)):
            if (
                face2.Surface.Axis.isEqual(face.Surface.Axis, 1e-5)
                and face2.Surface.Radius == rad
                and is_in_line(face2.Surface.Center, face.Surface.Axis, face.Surface.Center)
            ):
                face_index_0.append(i)

    # prueba same_faces, parece ok
    Faces_p = []
    for ind in face_index_0:
        Faces_p.append(solid.Faces[ind])

    face_index = [face_index_0[0]]  # la face de entrada

    for k in same_faces(Faces_p, tolerances):
        face_index.append(face_index_0[k])

    # commented to let plane cut rounded corner
    # if len(face_index_0)==len(face_index): #Esto evita un corte creo!!! no debería ser así en conversion
    #    return None

    for ind in reversed(face_index[1:]):
        if solid.Faces[ind].Area <= 1e-3:
            face_index.remove(ind)

        else:
            face2 = solid.Faces[ind]
            try:
                face2.tessellate(0.1)
                UVNodes.append(face2.getUVNodes())
            except RuntimeError:
                PR = face2.ParameterRange
                UVNode1 = (PR[0], PR[2])
                UVNode2 = (PR[1], PR[3])
                UVNodes.append([UVNode1, UVNode2])

    AngleRange = 0.0

    Uval = []
    for index in face_index:
        Range = solid.Faces[index].ParameterRange
        AngleRange = AngleRange + abs(Range[1] - Range[0])
        # if not(Range[0] in Uval) and not(Range[1] in Uval):
        Uval.append(Range[0])
        Uval.append(Range[1])

    if twoPi - AngleRange < 1.0e-2 or AngleRange < 1.0e-2:
        return None

    Uval_str_cl = []

    for i, elem1 in enumerate(Uval):
        num_str1 = f"{elem1:11.4E}"
        if abs(elem1) < 1.0e-5:
            num_str1 = "%11.4E" % 0.0

        if not (is_duplicate_in_list(num_str1, i, Uval)):
            Uval_str_cl.append(num_str1)

    if len(Uval_str_cl) < 2:
        logger.info("gen_plane_cylinder : Uval_str_cl should no be void")
        return None

    face_index_2 = [face_index[0], face_index[0]]

    Node_min = UVNodes[0][0]
    Node_max = UVNodes[0][1]

    dif1_0 = abs(float(Uval_str_cl[0]) - Node_min[0])
    dif2_0 = abs(float(Uval_str_cl[1]) - Node_max[0])

    # searching for minimum and maximum angle points

    for j, Nodes in enumerate(UVNodes):
        for elem in Nodes:
            dif1 = abs(float(Uval_str_cl[0]) - elem[0])
            dif2 = abs(float(Uval_str_cl[1]) - elem[0])

            if dif1 < dif1_0:
                Node_min = elem
                face_index_2[0] = face_index[j]
                dif1_0 = dif1
            if dif2 < dif2_0:
                Node_max = elem
                face_index_2[1] = face_index[j]
                dif2_0 = dif2

    V1 = solid.Faces[face_index_2[0]].valueAt(Node_min[0], Node_min[1])
    V2 = solid.Faces[face_index_2[1]].valueAt(Node_max[0], Node_max[1])

    if V1.isEqual(V2, 1e-5):
        logger.error("in the additional plane definition")
        return None

    normal = V2.sub(V1).cross(face.Surface.Axis)

    plane = Part.Plane(V1, normal).toShape()

    return plane


# Tolerance in this function are not the general once
# function should be reviewed
def gen_plane_cone(face, solid, tolerances):

    if face.Area < 1e-2:
        return None

    UVNodes = []
    face_index_0 = [solid.Faces.index(face)]

    try:
        face.tessellate(0.1)
        UVNodes.append(face.getUVNodes())
    except RuntimeError:
        PR = face.ParameterRange
        UVNode1 = (PR[0], PR[2])
        UVNode2 = (PR[1], PR[3])
        UVNodes.append([UVNode1, UVNode2])

    for i, face2 in enumerate(solid.Faces):

        if str(face2.Surface) == "<Cone object>" and not (face2.isEqual(face)):

            if (
                face2.Surface.Axis.isEqual(face.Surface.Axis, 1e-5)
                and face2.Surface.Apex.isEqual(face.Surface.Apex, 1e-5)
                and (face2.Surface.SemiAngle - face.Surface.SemiAngle) < 1e-6
            ):

                face_index_0.append(i)

    Faces_p = []
    for ind in face_index_0:
        Faces_p.append(solid.Faces[ind])

    face_index = [face_index_0[0]]  # la face de entrada

    for k in same_faces(Faces_p, tolerances):
        face_index.append(face_index_0[k])

    # same as cylinder commennt
    # if len(face_index_0)==len(face_index):
    #    return None

    for ind in reversed(face_index[1:]):
        if solid.Faces[ind].Area <= 1e-3:
            face_index.remove(ind)
        else:
            face2 = solid.Faces[ind]
            try:
                face2.tessellate(0.1)
                UVNodes.append(face2.getUVNodes())
            except RuntimeError:
                PR = face2.ParameterRange
                UVNode1 = (PR[0], PR[2])
                UVNode2 = (PR[1], PR[3])
                UVNodes.append([UVNode1, UVNode2])

    AngleRange = 0.0

    Uval = []

    for index in face_index:
        Range = solid.Faces[index].ParameterRange
        AngleRange = AngleRange + abs(Range[1] - Range[0])
        Uval.append(Range[0])
        Uval.append(Range[1])
    if twoPi - AngleRange < 1.0e-2 or AngleRange < 1e-2:
        return None

    Uval_str_cl = []

    for i, elem1 in enumerate(Uval):
        num_str1 = f"{elem1:11.4E}"
        if abs(elem1) < 1.0e-5:
            num_str1 = "%11.4E" % 0.0
        if not (is_duplicate_in_list(num_str1, i, Uval)):
            Uval_str_cl.append(num_str1)

    if len(Uval_str_cl) < 2:
        logger.info("gen_plane_cone : Uval_str_cl should no be void")
        return None

    face_index_2 = [face_index[0], face_index[0]]

    Node_min = UVNodes[0][0]
    Node_max = UVNodes[0][1]
    dif1_0 = abs(float(Uval_str_cl[0]) - Node_min[0])
    dif2_0 = abs(float(Uval_str_cl[1]) - Node_max[0])

    # searching for minimum and maximum angle points
    for j, Nodes in enumerate(UVNodes):
        for elem in Nodes:
            dif1 = abs(float(Uval_str_cl[0]) - elem[0])
            dif2 = abs(float(Uval_str_cl[1]) - elem[0])

            if dif1 < dif1_0:
                Node_min = elem
                face_index_2[0] = face_index[j]
                dif1_0 = dif1
            if dif2 < dif2_0:
                Node_max = elem
                face_index_2[1] = face_index[j]
                dif2_0 = dif2

    V1 = solid.Faces[face_index_2[0]].valueAt(Node_min[0], Node_min[1])
    V2 = solid.Faces[face_index_2[1]].valueAt(Node_max[0], Node_max[1])

    if V1.isEqual(V2, 1e-5):
        logger.error("in the additional plane definition")
        return None

    # normal=V2.sub(V1).cross(face.Surface.Axis)

    # plane=Part.Plane(V1,normal).toShape()
    plane = Part.Plane(V1, V2, face.Surface.Apex).toShape()

    return plane


def plane_2nd_order(solid_GU, face, flag_inv, tolerances, convex=True):
    planes = []

    if face is None:
        for face in solid_GU.Faces:
            if flag_inv:
                orient_temp = face.Orientation
                if orient_temp == "Forward":
                    orient = "Reversed"
                elif orient_temp == "Reversed":
                    orient = "Forward"
            else:
                orient = face.Orientation

            if convex and orient == "Reversed":  # convex face condition
                pp = None
                if str(face.Surface) == "<Cylinder object>":
                    pp = gen_plane_cylinder(face, solid_GU, tolerances)
                elif str(face.Surface) == "<Cone object>":
                    pp = gen_plane_cone(face, solid_GU, tolerances)
                if pp is not None:
                    planes.append(pp)
            # elif str(face.Surface)[0:6] == 'Sphere':
            #    return gen_plane_sphere(face,solid_GU)

    else:
        if flag_inv:
            orient_temp = face.Orientation
            if orient_temp == "Forward":
                orient = "Reversed"
            elif orient_temp == "Reversed":
                orient = "Forward"
        else:
            orient = face.Orientation

        if not convex and orient == "Forward":
            if str(face.Surface) == "<Cylinder object>":
                pp = gen_plane_cylinder(face, solid_GU, tolerances)
            elif str(face.Surface) == "<Cone object>":
                pp = gen_plane_cone(face, solid_GU, tolerances)
            if pp is not None:
                planes.append(pp)

    return planes


def split_planes(Solids, universe_box, options, tolerances, numeric_format):
    Bases = Solids[:]
    simpleSolid = []

    # Cut with real planes defining solid faces
    while True:
        newBases = []
        for base in Bases:
            cut_solids = split_p_planes_new(base, universe_box, options, tolerances, numeric_format)
            if len(cut_solids) == 1:
                simpleSolid.extend(cut_solids)
            else:
                newBases.extend(cut_solids)
        if len(newBases) == 0:
            break
        else:
            Bases = newBases

    return simpleSolid, 0


class ParallelPlanes:
    def __init__(self, plane):
        self.Axis = plane.Surf.Axis
        self.elements = [plane]
        self.count = 1

    def append(self, plane):
        if is_parallel(plane.Surf.Axis, self.Axis, 0.1):
            self.elements.append(plane)
            self.count += 1
            return True
        else:
            return False

    def sort(self):
        pos = []
        for i, p in enumerate(self.elements):
            d = self.Axis.dot(p.Surf.Position)
            pos.append((d, i))
        pos.sort()

        sort = []
        for p in pos:
            sort.append(self.elements[p[1]])
        self.elements = sort


def sort_planes(PlaneList, sortElements=False):
    if not PlaneList:
        return []
    newList = [ParallelPlanes(PlaneList[0])]
    for p in PlaneList[1:]:
        found = False
        for pp in newList:
            if pp.append(p):
                found = True
                break
        if not found:
            newList.append(ParallelPlanes(p))

    lenList = []
    for i, pp in enumerate(newList):
        if sortElements:
            pp.sort()
        lenList.append((pp.count, i))
    lenList.sort()

    sortedPlanes = []
    for l in lenList:
        sortedPlanes.append(newList[l[1]].elements)

    return sortedPlanes


def split_p_planes_new(solid, universe_box, options, tolerances, numeric_format):
    SPlanes = extract_surfaces(solid, "Planes", universe_box, options, tolerances, numeric_format, False)

    Planes = []
    for P in ("PX", "PY", "PZ", "P"):
        Planes.extend(SPlanes[P])

    Planes = sort_planes(Planes)

    if len(Planes) == 0:
        return [solid]
    if len(Planes[-1]) < options.nPlaneReverse:
        Planes.reverse()
    out_solid = [solid]
    for pp in Planes:
        for p in pp:
            p.build_surface()
        tools = tuple(p.shape for p in pp)
        comsolid = UF.split_bop(solid, tools, options.splitTolerance, options)

        if len(comsolid.Solids) > 1:
            out_solid = comsolid.Solids
            break
    return out_solid


def split_2nd_order(Solids, universe_box, options, tolerances, numeric_format):
    err = 0
    Base = Solids
    for kind in ["Cyl", "Cone", "Sph", "Tor"]:
        kindBase = []
        while True:
            cutBase = []
            for solid in Base:
                Surfaces = extract_surfaces(
                    solid,
                    kind,
                    universe_box,
                    options,
                    tolerances,
                    numeric_format,
                    False,
                )
                if len(Surfaces[kind]) == 0:
                    kindBase.append(solid)
                else:
                    simple = True
                    for s in Surfaces[kind]:
                        s.build_surface()
                        try:
                            comsolid = UF.split_bop(solid, [s.shape], options.splitTolerance, options)
                            solidsInCom = []
                            for s in comsolid.Solids:
                                if s.Volume > 1e-9:
                                    solidsInCom.append(s)

                            if len(solidsInCom) > 1:
                                simple = False
                                cutBase.extend(solidsInCom)
                                break
                        except:
                            logger.info("Failed split base with {s.shape.Faces[0].Surface} surface")
                            err += 2

                    if simple:
                        kindBase.append(solid)

            if len(cutBase) == 0:
                Base = kindBase
                break
            else:
                Base = cutBase

    return Base, err


def split_2nd_order_planes(Solids, options, tolerances):
    err = 0
    simpleSolid = []
    Bases = Solids
    while True:
        newBases = []
        for base in Bases:
            cut_solids, err = split_2nd_o_plane(base, options, tolerances)
            if len(cut_solids) == 1:
                simpleSolid.extend(cut_solids)
            else:
                newBases.extend(cut_solids)
        if len(newBases) == 0:
            break
        else:
            Bases = newBases

    return simpleSolid, err


def split_2nd_o_plane(solid, options, tolerances):

    err = 0
    flag_inv = CD.is_inverted(solid)
    solid_GU = GU.SolidGu(solid, tolerances=tolerances)
    planes = plane_2nd_order(solid_GU, None, flag_inv, tolerances, convex=True)
    if not planes:
        return [solid], err

    for p in planes:
        comsolid = UF.split_bop(solid, [p], options.splitTolerance, options)
        if not comsolid.Solids:
            comsolid = solid
            continue
        if len(comsolid.Solids) > 1:
            removed, err = remove_solids(comsolid.Solids)
            comsolid = Part.makeCompound(removed)
            if len(removed) > 1:
                break

    return comsolid.Solids, err


def remove_solids(Solids):
    err = 0
    Vol_tol = 1e-2
    Solids_Clean = []
    for solid in Solids:
        if solid.Volume <= Vol_tol:
            logger.warning(f"remove_solids degenerated solids are produced {solid.Volume}")
            err = 2
            continue
        Solids_Clean.append(solid)

    return Solids_Clean, err


def split_component(solidShape, universe_box, options, tolerances, numeric_format):
    err = 0
    err2 = 0

    Volini = solidShape.Volume
    Solids = solidShape.Solids
    # Split with explicit planes bounding the solid and
    # implicit planes interface of two 2nd order surfaces
    split0, err = split_planes(Solids, universe_box, options, tolerances, numeric_format)
    # Split with explicit 2nd order surfaces bounding the solid

    split1, err1 = split_2nd_order(split0, universe_box, options, tolerances, numeric_format)
    err += err1

    split, err2 = split_2nd_order_planes(split1, options, tolerances)
    err += err2
    Pieces = []
    for part in split:
        if part.Volume <= 1e-10:
            logger.warning(
                f"split_component degenerated solids are produced {part.Volume}",
            )
            err += 2
            continue
        Pieces.append(part)

    comsolid = Part.makeCompound(Pieces)

    Volend = comsolid.Volume
    Volch = (Volini - Volend) / Volini

    if abs(Volch) > 1e-2:  # 1% of Volume change
        logger.warning(f"Volume has changed after decomposition: {Volch:11.4E}")
        err += 4

    return comsolid, err


# TODO rename this function as there are two with the name name
def SplitSolid(solidShape, universe_box, options, tolerances, numeric_format):

    solid_parts = []

    for solid in solidShape.Solids:

        explode = split_full_cylinder(solid, options, tolerances, numeric_format)
        piece, err = split_component(explode, universe_box, options, tolerances, numeric_format)
        solid_parts.append(piece)

    return Part.makeCompound(solid_parts), err
