############################
# Module for Cell definiton #
#############################
import math

import FreeCAD
import Part

from ..Utils import BasicFunctions_part2 as BF
from ..Utils import Functions as UF
from ..Utils import Geometry_GU as GU
from ..Utils.BasicFunctions_part1 import (
    is_in_line,
    is_opposite,
    is_parallel,
    is_same_value,
    sign_plane,
)
from ..Utils.booleanFunction import BoolSequence, insert_in_sequence
from ..Utils.BooleanSolids import build_c_table_from_solids, remove_extra_surfaces
from ..Utils.Functions import GeounedSurface


def get_id(facein, surfaces, tolerances, options, numeric_format):

    surfin = str(facein)
    if surfin == "<Plane object>":
        if is_parallel(facein.Axis, FreeCAD.Vector(1, 0, 0), tolerances.pln_angle):
            p = "PX"
        elif is_parallel(facein.Axis, FreeCAD.Vector(0, 1, 0), tolerances.pln_angle):
            p = "PY"
        elif is_parallel(facein.Axis, FreeCAD.Vector(0, 0, 1), tolerances.pln_angle):
            p = "PZ"
        else:
            p = "P"

        for s in surfaces[p]:
            if BF.is_same_plane(
                facein,
                s.Surf,
                options=options,
                tolerances=tolerances,
                numeric_format=numeric_format,
                dtol=tolerances.pln_distance,
                atol=tolerances.pln_angle,
                rel_tol=tolerances.relativeTol,
            ):
                return s.Index

    elif surfin == "<Cylinder object>":
        for s in surfaces["Cyl"]:
            if BF.is_same_cylinder(
                facein,
                s.Surf,
                tolerances=tolerances,
                options=options,
                numeric_format=numeric_format,
                dtol=tolerances.cyl_distance,
                atol=tolerances.cyl_angle,
                rel_tol=tolerances.relativeTol,
            ):
                return s.Index

    elif surfin == "<Cone object>":
        for s in surfaces["Cone"]:
            if BF.is_same_cone(
                facein,
                s.Surf,
                dtol=tolerances.kne_distance,
                atol=tolerances.kne_angle,
                rel_tol=tolerances.relativeTol,
            ):
                return s.Index

    elif surfin[0:6] == "Sphere":
        for s in surfaces["Sph"]:
            if BF.is_same_sphere(
                sph1=facein,
                sph2=s.Surf,
                tolerance=tolerances.sph_distance,
                rel_tol=tolerances.relativeTol,
            ):
                return s.Index

    elif surfin == "<Toroid object>":
        for s in surfaces["Tor"]:
            if BF.is_same_torus(
                facein,
                s.Surf,
                dtol=tolerances.tor_distance,
                atol=tolerances.tor_angle,
                rel_tol=tolerances.relativeTol,
            ):
                return s.Index

    return 0


def is_inverted(solid):

    face = solid.Faces[0]

    # u=(face.Surface.bounds()[0]+face.Surface.bounds()[1])/2.0 # entre 0 y 2pi si es completo
    # v=face.Surface.bounds()[0]+(face.Surface.bounds()[3]-face.Surface.bounds()[2])/3.0 # a lo largo del eje
    parameter_range = face.ParameterRange
    u = (parameter_range[1] + parameter_range[0]) / 2.0
    v = (parameter_range[3] + parameter_range[2]) / 2.0

    if str(face.Surface) == "<Cylinder object>":
        dist1 = face.Surface.value(u, v).distanceToLine(
            face.Surface.Center, face.Surface.Axis
        )
        dist2 = (
            face.Surface.value(u, v)
            .add(face.Surface.normal(u, v).multiply(1.0e-6))
            .distanceToLine(face.Surface.Center, face.Surface.Axis)
        )
        if (dist2 - dist1) < 0.0:
            # The normal of the cylinder is going inside
            return True
    elif str(face.Surface) == "<Cone object>":
        dist1 = face.Surface.value(u, v).distanceToLine(
            face.Surface.Apex, face.Surface.Axis
        )
        dist2 = (
            face.Surface.value(u, v)
            .add(face.Surface.normal(u, v).multiply(1.0e-6))
            .distanceToLine(face.Surface.Apex, face.Surface.Axis)
        )
        if (dist2 - dist1) < 0.0:
            # The normal of the cylinder is going inside
            return True
    # MIO
    elif str(face.Surface)[0:6] == "Sphere":
        # radii = point - center
        radii = face.Surface.value(u, v).add(face.Surface.Center.multiply(-1))
        radii_b = (
            face.Surface.value(u, v)
            .add(face.Surface.normal(u, v).multiply(1.0e-6))
            .add(face.Surface.Center.multiply(-1))
        )
        # radii_b  = radii.add( face.Surface.normal(u,v).multiply(1.0e-6) )
        if (radii_b.Length - radii.Length) < 0.0:
            # An increasing of the radii vector in the normal direction decreases the radii: oposite normal direction
            return True

    elif str(face.Surface) == "<Plane object>":
        dist1 = face.CenterOfMass.distanceToPoint(solid.BoundBox.Center)
        dist2 = face.CenterOfMass.add(
            face.normalAt(u, v).multiply(1.0e-6)
        ).distanceToPoint(solid.BoundBox.Center)
        point2 = face.CenterOfMass.add(face.normalAt(u, v).multiply(1.0e-6))
        if solid.isInside(point2, 1e-7, False):
            return True

    return False


def gen_plane(face, solid, verbose):
    """Generate an additional plane when convex surfaces of second order are presented as a face of the solid"""

    surf = face.Surface
    if str(surf) == "<Cylinder object>":
        return gen_plane_cylinder(face, solid, verbose)
    if str(surf) == "<Cone object>":
        return gen_plane_cone(face, solid)
    if str(surf) == "Sphere":
        return gen_plane_sphere(face, solid)


def get_closed_ranges(solid, face_index):

    u_nodes = []
    for index in face_index:
        URange = solid.Faces[index].ParameterRange
        u_nodes.append((URange[0], index))
        u_nodes.append((URange[1], index))
    u_nodes.sort()

    closed_range = get_intervals(u_nodes)

    a_min = closed_range[0][0][0]
    a_max = closed_range[-1][1][0]

    if abs(a_max - a_min - 2.0 * math.pi) < 1e-2:
        if len(closed_range) == 1:
            closed_face = True
        else:
            endPoint = (closed_range[-1][0][0] - 2 * math.pi, closed_range[-1][0][1])
            closed_range[0][0] = endPoint
            closed_range[0][2].update(closed_range[-1][2])
            del closed_range[-1]

            if len(closed_range) == 1:
                if (
                    abs(closed_range[0][1][0] - closed_range[0][0][0] - 2.0 * math.pi)
                    < 1e-2
                ):
                    closed_face = True
                else:
                    closed_face = False
            else:
                closed_face = False
    else:
        closed_face = False
    return closed_range, closed_face


def get_intervals(u_nodes):
    closed_ranges = []
    pos_min = dict()
    pos_max = dict()
    for i, node in enumerate(u_nodes):
        if node[1] not in pos_min.keys():
            pos_min[node[1]] = i
        else:
            pos_max[node[1]] = i

    u_min = u_nodes[0]
    i_pos = pos_max[u_min[1]]

    while True:
        x = u_nodes[i_pos]
        end = True
        for i in range(i_pos + 1, len(u_nodes)):
            mxt_int = u_nodes[i][1]
            if (
                u_nodes[pos_min[mxt_int]][0] - x[0]
            ) < 1e-5:  # x pos is > min boundary of the next inteval inside precision 1e-5
                i_pos = pos_max[mxt_int]
                end = False
                break

        if end:
            u_max = x
            closed_ranges.append([u_min, u_max])
            i_pos += 1
            if i_pos < len(u_nodes):
                u_min = u_nodes[i_pos]
                i_pos = pos_max[u_min[1]]
            else:
                break

    for closed_range in closed_ranges:
        index = set()
        xmin = closed_range[0][0]
        xmax = closed_range[1][0]
        for interval in u_nodes:
            x = interval[0]
            if (xmin - x) < 1.0e-5 and (x - xmax) < 1.0e-5:
                index.add(interval[1])
        closed_range.append(index)

    return closed_ranges


def get_u_value_boundary(solid, face_index, my_index):

    face_u_ranges, closed_face = get_closed_ranges(solid, face_index)
    if closed_face:
        return None, None

    for face_u_range in face_u_ranges:
        if my_index in face_u_range[2]:
            u_min, u_max = face_u_range[0:2]
            return u_min, u_max


def gen_plane_sphere(face, solid):
    same_faces = []
    same_faces.append(face)

    for f in solid.Faces:
        if f.isEqual(face) or str(f.Surface) != "Sphere":
            continue
        if (
            f.Surface.Center == face.Surface.Center
            and f.Surface.Radius == face.Surface.Radius
        ):
            # print 'Warning: coincident sphere faces are the same'
            for f2 in same_faces:
                if f.distToShape(f2)[0] < 1e-6:
                    same_faces.append(f)
                    break

    # print same_faces
    normal = FreeCAD.Vector(0, 0, 0)
    for f in same_faces:
        normal += f.Area * (f.CenterOfMass - face.Surface.Center)

    return Part.Plane(face.Surface.Center, normal).toShape()


def gen_plane_cylinder(face, solid, verbose):

    surf = face.Surface
    rad = surf.Radius

    if str(surf) != "<Cylinder object>":
        return None

    my_index = solid.Faces.index(face)
    face_index = [my_index]

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < min_area:
            if verbose:
                print(
                    f"Warning: {str(surf)} surface removed from cell definition. Face area < Min area ({face2.Area} < {min_area}) "
                )
            continue
        if str(face2.Surface) == "<Cylinder object>" and not (face2.isEqual(face)):
            if (
                face2.Surface.Axis.isEqual(face.Surface.Axis, 1e-5)
                and face2.Surface.Radius == rad
                and is_in_line(
                    face2.Surface.Center, face.Surface.Axis, face.Surface.Center
                )
            ):
                # print 'Warning: coincident cylinder faces are the same'
                face_index.append(i)

    u_min, u_max = get_u_value_boundary(solid, face_index, my_index)
    if u_min is None:
        return None

    u_1, i1 = u_min
    u_2, i2 = u_max

    v_1 = solid.Faces[i1].ParameterRange[2]
    v_2 = solid.Faces[i2].ParameterRange[2]

    p1 = solid.Faces[i1].valueAt(u_1, v_1)
    p2 = solid.Faces[i2].valueAt(u_2, v_2)

    if p1.isEqual(p2, 1e-5):
        if verbose:
            print("Error in the additional place definition")
        return None

    normal = p2.sub(p1).cross(face.Surface.Axis)
    plane = Part.Plane(p1, normal).toShape()

    return plane


def gen_plane_cylinder_old(face, solid, verbose):

    surf = face.Surface
    rad = surf.Radius

    if str(surf) != "<Cylinder object>":
        return None

    face_index = [solid.Faces.index(face)]

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < min_area:
            if verbose:
                print(
                    f"Warning: {str(surf)} surface removed from cell definition. Face area < Min area ({face2.Area} < {min_area}) "
                )
            continue
        if str(face2.Surface) == "<Cylinder object>" and not (face2.isEqual(face)):
            if (
                face2.Surface.Axis.isEqual(face.Surface.Axis, 1e-5)
                and face2.Surface.Radius == rad
                and is_in_line(
                    face2.Surface.Center, face.Surface.Axis, face.Surface.Center
                )
            ):
                # print 'Warning: coincident cylinder faces are the same'
                face_index.append(i)

    angle_range = 0.0
    U_val = []
    for index in face_index:
        Range = solid.Faces[index].ParameterRange
        angle_range = angle_range + abs(Range[1] - Range[0])
        if not (Range[0] in U_val) and not (Range[1] in U_val):
            U_val.append(Range[0])
            U_val.append(Range[1])
    if 2.0 * math.pi - angle_range < 1e-2:
        return None

    UVNodes = []
    for index in face_index:
        face2 = solid.Faces[index]
        try:
            UVNodes.append(face2.getUVNodes())
        except RuntimeError:
            UVNodes.append(face2.getUVNodes())

    U_val_str_cl = []
    for i, elem1 in enumerate(U_val):
        num_str1 = f"{elem1:11.4E}"
        if abs(elem1) < 1.0e-5:
            num_str1 = "%11.4E" % 0.0
        if not (BF.is_duplicate_in_list(num_str1, i, U_val)):
            U_val_str_cl.append(num_str1)

    face_index_2 = [face_index[0], face_index[0]]

    node_min = UVNodes[0][0]
    node_max = UVNodes[0][1]

    dif1_0 = abs(float(U_val_str_cl[0]) - node_min[0])
    dif2_0 = abs(float(U_val_str_cl[1]) - node_max[0])

    # searching for minimum and maximum angle points
    for j, Nodes in enumerate(UVNodes):
        for elem in Nodes:
            dif1 = abs(float(U_val_str_cl[0]) - elem[0])
            dif2 = abs(float(U_val_str_cl[1]) - elem[0])

            if dif1 < dif1_0:
                node_min = elem
                face_index_2[0] = face_index[j]
                dif1_0 = dif1
            if dif2 < dif2_0:
                node_max = elem
                face_index_2[1] = face_index[j]
                dif2_0 = dif2

    v_1 = solid.Faces[face_index_2[0]].valueAt(node_min[0], node_min[1])
    v_2 = solid.Faces[face_index_2[1]].valueAt(node_max[0], node_max[1])

    if v_1.isEqual(v_2, 1e-5):
        if verbose:
            print("Error in the additional place definition")
        return None

    normal = v_2.sub(v_1).cross(face.Surface.Axis)
    plane = Part.Plane(v_1, normal).toShape()

    return plane


def gen_plane_cone(face, solid):

    Surf = face.Surface
    if str(Surf) != "<Cone object>":
        return None

    myIndex = solid.Faces.index(face)
    face_index = [myIndex]

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < min_area:
            if verbose:
                print(
                    f"Warning: {str(Surf)} surface removed from cell definition. Face area < Min area ({face2.Area} < {min_area}) "
                )
            continue
        if str(face2.Surface) == "<Cone object>" and not (face2.isEqual(face)):
            if (
                face2.Surface.Axis.isEqual(face.Surface.Axis, 1e-5)
                and face2.Surface.Apex.isEqual(face.Surface.Apex, 1e-5)
                and (face2.Surface.SemiAngle - face.Surface.SemiAngle) < 1e-6
            ):
                face_index.append(i)

    u_min, u_max = get_u_value_boundary(solid, face_index, myIndex)
    if u_min is None:
        return None

    u_1, i1 = u_min
    u_2, i2 = u_max

    v_1 = solid.Faces[i1].ParameterRange[2]
    v_2 = solid.Faces[i2].ParameterRange[2]

    p1 = solid.Faces[i1].valueAt(u_1, v_1)
    p2 = solid.Faces[i2].valueAt(u_2, v_2)

    if p1.isEqual(p2, 1e-5):
        if verbose:
            print("Error in the additional place definition")
        return None

    plane = Part.Plane(p1, p2, face.Surface.Apex).toShape()

    return plane


def gen_plane_cone_old(face, solid, verbose):

    surf = face.Surface
    if str(surf) != "<Cone object>":
        return None

    face_index = [solid.Faces.index(face)]

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < min_area:
            if verbose:
                print(
                    f"Warning: {str(surf)} surface removed from cell definition. Face area < Min area ({face2.Area} < {min_area}) "
                )
            continue
        if str(face2.Surface) == "<Cone object>" and not (face2.isEqual(face)):
            if (
                face2.Surface.Axis.isEqual(face.Surface.Axis, 1e-5)
                and face2.Surface.Apex.isEqual(face.Surface.Apex, 1e-5)
                and (face2.Surface.SemiAngle - face.Surface.SemiAngle) < 1e-6
            ):
                face_index.append(i)

    angle_range = 0.0
    u_val = []
    for index in face_index:
        parameter_range = solid.Faces[index].ParameterRange
        angle_range = angle_range + abs(parameter_range[1] - parameter_range[0])
        u_val.append(parameter_range[0])
        u_val.append(parameter_range[1])

    if 2.0 * math.pi - angle_range < 1e-2:
        return None

    uv_nodes = []
    for index in face_index:
        face2 = solid.Faces[index]
        try:
            uv_nodes.append(face2.getUVNodes())
        except RuntimeError:
            face.tessellate(1.0, True)
            uv_nodes.append(face2.getUVNodes())

    u_val_str_cl = []

    for i, elem1 in enumerate(u_val):
        num_str1 = f"{elem1:11.4E}"
        if abs(elem1) < 1.0e-5:
            num_str1 = "%11.4E" % 0.0
        if not (BF.is_duplicate_in_list(num_str1, i, u_val)):
            u_val_str_cl.append(num_str1)

    face_index_2 = [face_index[0], face_index[0]]

    node_min = uv_nodes[0][0]
    node_max = uv_nodes[0][1]
    dif1_0 = abs(float(u_val_str_cl[0]) - node_min[0])
    dif2_0 = abs(float(u_val_str_cl[1]) - node_max[0])

    # searching for minimum and maximum angle points
    for j, Nodes in enumerate(uv_nodes):
        for elem in Nodes:
            dif1 = abs(float(u_val_str_cl[0]) - elem[0])
            dif2 = abs(float(u_val_str_cl[1]) - elem[0])
            if dif1 < dif1_0:
                node_min = elem
                face_index_2[0] = face_index[j]
                dif1_0 = dif1
            if dif2 < dif2_0:
                node_max = elem
                face_index_2[1] = face_index[j]
                dif2_0 = dif2

    v_1 = solid.Faces[face_index_2[0]].valueAt(node_min[0], node_min[1])
    v_2 = solid.Faces[face_index_2[1]].valueAt(node_max[0], node_max[1])

    if v_1.isEqual(v_2, 1e-5):
        if verbose:
            print("Error in the additional place definition")
        return None

    plane = Part.Plane(v_1, v_2, face.Surface.Apex).toShape()

    return plane


def gen_torus_annex_u_planes(face, u_params, tolerances):

    if is_parallel(face.Surface.Axis, FreeCAD.Vector(1, 0, 0), tolerances.tor_angle):
        axis = FreeCAD.Vector(1, 0, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 1, 0), tolerances.tor_angle):
        axis = FreeCAD.Vector(0, 1, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 0, 1), tolerances.tor_angle):
        axis = FreeCAD.Vector(0, 0, 1)

    center = face.Surface.Center
    p1 = face.valueAt(u_params[0], 0.0)
    p2 = face.valueAt(u_params[1], 0.0)
    pmid = face.valueAt(0.5 * (u_params[0] + u_params[1]), 0.0)

    if is_same_value(abs(u_params[1] - u_params[0]), math.pi, tolerances.value):
        d = axis.cross(p2 - p1)
        d.normalize()
        if d.dot(pmid - center) < 0:
            d = -d
        return (
            (center, d, face.Surface.MajorRadius, face.Surface.MajorRadius),
            None,
        ), False

    elif u_params[1] - u_params[0] < math.pi:
        d = axis.cross(p2 - p1)
        d.normalize()
        if d.dot(pmid - center) < 0:
            d = -d
        return (
            (center, d, face.Surface.MajorRadius, face.Surface.MajorRadius),
            None,
        ), False

    else:
        d1 = axis.cross(p1)
        d1.normalize()
        if d1.dot(pmid - center) < 0:
            d1 = -d1

        d2 = axis.cross(p2)
        d2.normalize()
        if d2.dot(pmid - center) < 0:
            d2 = -d2

        return (
            (center, d1, face.Surface.MajorRadius, face.Surface.MajorRadius),
            (center, d2, face.Surface.MajorRadius, face.Surface.MajorRadius),
        ), True  # (d1 : d2)


def gen_torus_annex_u_planes_org(face, u_params):

    if is_parallel(face.Surface.Axis, FreeCAD.Vector(1, 0, 0), tor_angle):
        axis = FreeCAD.Vector(1, 0, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 1, 0), tor_angle):
        axis = FreeCAD.Vector(0, 1, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 0, 1), tor_angle):
        axis = FreeCAD.Vector(0, 0, 1)

    center = face.Surface.Center
    p1 = face.valueAt(u_params[0], 0.0)
    p2 = face.valueAt(u_params[1], 0.0)
    pmid = face.valueAt(0.5 * (u_params[0] + u_params[1]), 0.0)

    if is_same_value(abs(u_params[1] - u_params[0]), math.pi, value):
        d = axis.cross(p2 - p1)
        d.normalize()
        if pmid.dot(d) < 0:
            d = -d
        return ((center, d, face.Surface.MajorRadius), None), False

    else:
        d1 = axis.cross(p1)
        d1.normalize()
        if pmid.dot(d1) < 0:
            d1 = -d1

        d2 = axis.cross(p2)
        d2.normalize()
        if pmid.dot(d2) < 0:
            d2 = -d2

        if u_params[1] - u_params[0] < math.pi:
            return (
                (center, d1, face.Surface.MajorRadius, face.Surface.MajorRadius),
                (center, d2, face.Surface.MajorRadius, face.Surface.MajorRadius),
            ), False  # ( d1 d2 )
        else:
            return (
                (center, d1, face.Surface.MajorRadius, face.Surface.MajorRadius),
                (center, d2, face.Surface.MajorRadius, face.Surface.MajorRadius),
            ), True  # (d1 : d2)


def gen_torus_annex_v_surface(face, v_params, tolerances, force_cylinder=False):
    if is_parallel(face.Surface.Axis, FreeCAD.Vector(1, 0, 0), tolerances.tor_angle):
        axis = FreeCAD.Vector(1, 0, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 1, 0), tolerances.tor_angle):
        axis = FreeCAD.Vector(0, 1, 0)
    elif is_parallel(face.Surface.Axis, FreeCAD.Vector(0, 0, 1), tolerances.tor_angle):
        axis = FreeCAD.Vector(0, 0, 1)

    p1 = face.valueAt(0.0, v_params[0]) - face.Surface.Center
    z1 = p1.dot(axis)
    d1 = p1.cross(axis).Length

    p2 = face.valueAt(0.0, v_params[1]) - face.Surface.Center
    z2 = p2.dot(axis)
    d2 = p2.cross(axis).Length

    if is_same_value(z1, z2, tolerances.distance):
        surf_type = "Plane"
        center = face.Surface.Center + z1 * axis
        v_mid = (v_params[0] + v_params[1]) * 0.5
        p_mid = face.valueAt(0, v_mid) - face.Surface.Center
        if p_mid.dot(axis) < z1:
            in_surf = True
        else:
            in_surf = False
        return (
            (center, axis, face.Surface.MajorRadius, face.Surface.MajorRadius),
            surf_type,
            in_surf,
        )

    elif is_same_value(d1, d2, tolerances.distance) or force_cylinder:
        surf_type = "Cylinder"
        radius = min(d1, d2)
        center = face.Surface.Center
        if is_same_value(d1, face.Surface.MajorRadius, tolerances.distance):
            v_mid = (v_params[0] + v_params[1]) * 0.5
            p_mid = face.valueAt(0, v_mid) - center
            if p_mid.cross(axis).Length < face.Surface.MajorRadius:
                in_surf = True
            v_mid = (v_params[0] + v_params[1]) * 0.5
            p_mid = face.valueAt(0, v_mid) - center
            if p_mid.cross(axis).Length < face.Surface.MajorRadius:
                in_surf = True
                radius = max(d1, d2)
            else:
                in_surf = False
        else:
            if d1 < face.Surface.MajorRadius:
                in_surf = True
                radius = max(d1, d2)
            else:
                in_surf = False
        return (center, axis, radius, face.Surface.MinorRadius), surf_type, in_surf

    else:
        surf_type = "Cone"
        za = (z2 * d1 - z1 * d2) / (d1 - d2)
        apex = face.Surface.Center + za * axis
        semi_angle = abs(math.atan(d1 / (z1 - za)))

        cone_axis = axis if (z1 - za) > 0.0 else -axis

        v_mid = (v_params[0] + v_params[1]) * 0.5
        p_mid = face.valueAt(0, v_mid) - face.Surface.Center
        z_mid = p_mid.dot(axis)
        d_mid = p_mid.cross(axis).Length

        d_cone = d1 * (z_mid - za) / (z1 - za)
        in_surf = True if d_mid < d_cone else False

        return (
            (
                apex,
                cone_axis,
                semi_angle,
                face.Surface.MinorRadius,
                face.Surface.MajorRadius,
            ),
            surf_type,
            in_surf,
        )


def cellDef(meta_obj, surfaces, universe_box, options, tolerances, numeric_format):

    solids = meta_obj.Solids
    del_list = []

    piece_def = BoolSequence(operator="OR")
    iece_obj = []
    cones = set()
    for isol, solid in enumerate(solids):
        surf_piece = []
        surf_obj = []
        extra_plane_reverse = dict()

        flag_inv = is_inverted(solid)
        solid_gu = GU.SolidGu(
            solid=solid,
            tolerances=tolerances,
            plane3Pts=False,
        )
        last_torus = -1
        for iface, face in enumerate(solid_gu.Faces):
            surface_type = str(face.Surface)
            if abs(face.Area) < tolerances.min_area:
                if options.verbose:
                    print(
                        f"Warning: {surface_type} surface removed from cell definition. Face area < Min area ({face.Area} < {min_area}) "
                    )
                continue
            if face.Area < 0:
                if options.verbose:
                    print("Warning : Negative surface Area")
            if face.Orientation not in ("Forward", "Reversed"):
                continue
            if flag_inv:
                orient_temp = face.Orientation
                if orient_temp == "Forward":
                    orient = "Reversed"
                elif orient_temp == "Reversed":
                    orient = "Forward"
            else:
                orient = face.Orientation

            if "Sphere" in surface_type:
                surface_type = "Sphere"

            # cone additional plane is added afterward
            if (
                surface_type in ("<Cylinder object>", "<Cone object>", "Sphere")
                and orient == "Reversed"
            ):
                # cone additional plane is added afterward
                id_face = get_id(
                    facein=face.Surface,
                    surfaces=surfaces,
                    tolerances=tolerances,
                    options=options,
                    numeric_format=numeric_format,
                )
                if surface_type == "<Cone object>":
                    cones.add(id_face)
                if str(id_face) not in surf_piece:
                    surf_piece.append(str(id_face))
                    surf_obj.append(face)

                try:
                    plane = gen_plane(face, solid_gu)
                    if plane is not None:
                        plane = GU.PlaneGu(plane)
                except:
                    plane = None
                    if options.verbose:
                        print("Warning: generation of additional plane has failed")

                if plane is not None:
                    p = GeounedSurface(
                        ("Plane", (plane.Position, plane.Axis, plane.dim1, plane.dim2)),
                        universe_box,
                        Face="Build",
                    )

                    id, exist = surfaces.add_plane(p)
                    sign = sign_plane(face.CenterOfMass, p)
                    if exist:
                        pp = surfaces.get_surface(id)
                        if is_opposite(p.Surf.Axis, pp.Surf.Axis, tolerances.angle):
                            id = -id
                    id *= sign

                    if id_face not in extra_plane_reverse.keys():
                        extra_plane_reverse[id_face] = [str(id)]
                        surf_obj.append(p.shape)
                    else:
                        if str(id) not in extra_plane_reverse[id_face]:
                            extra_plane_reverse[id_face].append(str(id))
                            surf_obj.append(p.shape)

            elif surface_type == "<Toroid object>":

                if (
                    is_parallel(
                        face.Surface.Axis, FreeCAD.Vector(1, 0, 0), tolerances.angle
                    )
                    or is_parallel(
                        face.Surface.Axis, FreeCAD.Vector(0, 1, 0), tolerances.angle
                    )
                    or is_parallel(
                        face.Surface.Axis, FreeCAD.Vector(0, 0, 1), tolerances.angle
                    )
                ):

                    idT = get_id(
                        facein=face.Surface,
                        surfaces=surfaces,
                        tolerances=tolerances,
                        options=options,
                        numeric_format=numeric_format,
                    )
                    index, u_params = solid_gu.TorusUParams[iface]
                    if index == last_torus:
                        continue
                    last_torus = index

                    # add if necesary additional planes following U variable
                    u_closed, u_minMax = u_params
                    # u_closed = True
                    if not u_closed:
                        planes, ORop = gen_torus_annex_u_planes(
                            face=face, u_params=u_minMax, tolerances=tolerances
                        )
                        plane1, plane2 = planes
                        plane = GeounedSurface(
                            ("Plane", plane1), universe_box, Face="Build"
                        )
                        id1, exist = surfaces.add_plane(
                            plane=plane,
                            tolerances=tolerances,
                            options=options,
                            numeric_format=numeric_format,
                        )
                        if exist:
                            p = surfaces.get_surface(id1)
                            if is_opposite(
                                vector_1=plane.Surf.Axis,
                                vector_2=p.Surf.Axis,
                                tolerance=tolerances.pln_angle,
                            ):
                                id1 = -id1

                        if plane2 is None:
                            u_var = "%i" % id1
                        else:
                            plane = GeounedSurface(
                                ("Plane", plane2), universe_box, Face="Build"
                            )
                            id2, exist = surfaces.add_plane(plane)
                            if exist:
                                p = surfaces.get_surface(id2)
                                if is_opposite(
                                    vector_1=plane.Surf.Axis,
                                    vector_2=p.Surf.Axis,
                                    tolerance=tolerances.pln_angle,
                                ):
                                    id2 = -id2

                            u_var = (
                                "(%i : %i)" % (id1, id2)
                                if ORop
                                else "%i %i" % (id1, id2)
                            )

                    else:
                        u_var = ""

                    # add if necesary additional surface following V variable
                    if orient == "Forward":
                        v_var = "-%i" % idT

                    else:
                        index, Vparams = solid_gu.TorusVParams[iface]
                        VClosed, VminMax = Vparams
                        if VClosed:
                            v_var = "%i" % idT
                        else:
                            surf_params, surf_type, in_surf = gen_torus_annex_v_surface(
                                face=face,
                                v_params=VminMax,
                                tolerances=tolerances,
                                force_cylinder=options.forceCylinder,
                            )

                            if surf_type == "Cone":
                                cone = GeounedSurface(
                                    ("Cone", surf_params), universe_box, Face="Build"
                                )
                                id2, exist = surfaces.add_cone(
                                    cone=cone, tolerances=tolerances
                                )

                            elif surf_type == "Cylinder":
                                cyl = GeounedSurface(
                                    ("Cylinder", surf_params),
                                    universe_box,
                                    Face="Build",
                                )
                                id2, exist = surfaces.add_cylinder(
                                    cyl=cyl,
                                    tolerances=tolerances,
                                    options=options,
                                    numeric_format=numeric_format,
                                )

                            elif surf_type == "Plane":
                                plane = GeounedSurface(
                                    ("Plane", surf_params), universe_box, Face="Build"
                                )
                                id2, exist = surfaces.add_plane(
                                    plane=plane,
                                    tolerances=tolerances,
                                    options=options,
                                    numeric_format=numeric_format,
                                )
                                if exist:
                                    p = surfaces.get_surface(id2)
                                    if is_opposite(
                                        plane.Surf.Axis,
                                        p.Surf.Axis,
                                        tolerances.pln_angle,
                                    ):
                                        id2 = -id2

                            v_var = "%i %i" % (idT, -id2 if in_surf else id2)

                    var = v_var if u_closed else " ".join((v_var, u_var))
                    if var not in surf_piece:
                        surf_piece.append(var)
                        surf_obj.append(face)
                else:
                    if options.verbose:
                        print(
                            "Only Torus with axis along X, Y , Z axis can be reproduced"
                        )
            else:
                id = get_id(
                    facein=face.Surface,
                    surfaces=surfaces,
                    tolerances=tolerances,
                    options=options,
                    numeric_format=numeric_format,
                )
                if surface_type == "<Cone object>":
                    cones.add(-id)

                surf = face
                if id == 0:
                    if options.verbose:
                        print("Warning: ", surface_type, " not found in surface list")
                    if surface_type == "<Plane object>":
                        dim1 = face.ParameterRange[1] - face.ParameterRange[0]
                        dim2 = face.ParameterRange[3] - face.ParameterRange[2]
                        plane = GeounedSurface(
                            (
                                "Plane",
                                (face.Surface.Position, face.Surface.Axis, dim1, dim2),
                            ),
                            universe_box,
                            Face="Build",
                        )
                        id, exist = surfaces.add_plane(plane)
                        surf = plane.shape
                    elif surface_type == "<Cylinder object>":
                        dim_l = face.ParameterRange[3] - face.ParameterRange[2]
                        cylinder = GeounedSurface(
                            (
                                "Cylinder",
                                (
                                    face.Surface.Center,
                                    face.Surface.Axis,
                                    face.Surface.Radius,
                                    dim_l,
                                ),
                            ),
                            universe_box,
                            Face="Build",
                        )
                        id, exist = surfaces.add_cylinder(cylinder)
                        surf = cylinder.shape

                if orient == "Reversed":
                    var = id
                elif orient == "Forward":
                    var = -id

                if surface_type == "<Plane object>":
                    s = surfaces.get_surface(id)
                    if is_opposite(
                        face.Surface.Axis, s.Surf.Axis, tolerances.pln_angle
                    ):
                        var = -var

                if str(var) in surf_piece:
                    continue

                surf_piece.append(str(var))
                surf_obj.append(surf)

        if extra_plane_reverse:
            for extra in extra_plane_reverse.values():
                if len(extra) == 1:
                    if extra[0] not in surf_piece:
                        surf_piece.append(extra[0])
                else:
                    surf_piece.append(f"({':'.join(extra)})")

        surf_piece_bool = BoolSequence(" ".join(surf_piece))
        # possible expresion for e
        #  i1
        #  i1 i2surf_piece
        #  i1 i2 i3
        #  (i1:i2)
        #  i1 (i2:i3)

        #         semi = e.find(':')
        #         blk  = e.find(' ')
        #         print (e)
        #         if semi != -1 :
        #           orTerm = expOR(int(e[1:semi]),int(e[semi+1:-1]))
        #           surf_piece_bool.add(orTerm)
        #         elif blk != -1 :
        #           surf_piece_bool.add(int(e[1:blk]),int(e[blk+1:-1]))
        #         else :
        #          surf_piece_bool.add(int(e))

        if surf_piece_bool.elements:
            piece_def.append(surf_piece_bool)
            iece_obj.append(surf_obj)
        else:
            del_list.append(isol)

    for isol in reversed(del_list):
        del meta_obj.Solids[isol]
    meta_obj.set_definition(piece_def)
    meta_obj.set_faces(iece_obj)
    return tuple(cones)


def get_surf_value(definition, reverse=False):

    if definition.level == 0:
        if reverse:
            surf = {-i for i in definition.elements}
        else:
            surf = set(definition.elements)
    else:
        surf = set()
        for e in definition.elements:
            if e.operator == "AND":
                if reverse:
                    surf = {-i for i in e.elements}
                else:
                    surf = set(e.elements)
                break
    return surf


def append_comp(new_cell, cell_def, cell_cad, meta_complementary):

    surf_cell = get_surf_value(cell_def, True)
    if meta_complementary.Definition.operator == "AND":
        if not cell_cad.BoundBox.intersect(meta_complementary.CADSolid.BoundBox):
            return False
        Seq = meta_complementary.Definition
        surfComp = get_surf_value(Seq, False)
        if len(surfComp & surf_cell) > 0:
            return False
        new_cell.append(Seq.get_complementary())
        return True

    else:
        append = False
        for i, compPart in enumerate(meta_complementary.Solids):
            if not cell_cad.BoundBox.intersect(compPart.BoundBox):
                continue
            Seq = meta_complementary.Definition.elements[i]
            surfComp = get_surf_value(Seq, False)
            if len(surfComp & surf_cell) > 0:
                continue
            append = True
            new_cell.append(Seq.get_complementary())
        return append


def no_overlapping_cell(metaList, surfaces, options):

    Surfs = {}
    for lst in surfaces.values():
        for s in lst:
            Surfs[s.Index] = s

    new_definition_list = []
    metaList[0].set_cad_solid()

    for i, m in enumerate(metaList[1:]):
        m.set_cad_solid()

        if m.Definition.operator == "AND":
            new_def = BoolSequence(operator="AND")
            new_def.append(m.Definition.copy())
            simplify = False
            for mm in metaList[: i + 1]:
                simp = append_comp(new_def, m.Definition, m.CADSolid, mm)
                if simp:
                    simplify = True
            simpTerm = [simplify]

        else:
            new_def = BoolSequence(operator="OR")
            simpTerm = []
            for j, partSolid in enumerate(m.Solids):
                subDef = BoolSequence(operator="AND")
                subDef.append(m.Definition.elements[j].copy())
                simplify = False
                for mm in metaList[: i + 1]:
                    simp = append_comp(subDef, m.Definition.elements[j], partSolid, mm)
                    if simp:
                        simplify = True
                simpTerm.append(simplify)
                new_def.append(subDef)
        new_definition_list.append((new_def, simpTerm))

    for m, t_def_and_simplify in zip(metaList[1:], new_definition_list):
        t_def, simplify = t_def_and_simplify
        if True in simplify:
            print(f"reduce cell {m.__id__}")
            box = UF.get_box(m)

            # evaluate only diagonal elements of the Constraint Table (fastest) and remove surface not
            # crossing in the solid boundBox
            CT = build_c_table_from_solids(
                Box=box,
                SurfInfo=(tuple(t_def.get_surfaces_numbers()), Surfs),
                options=options,
                option="diag",
            )

            new_def = remove_extra_surfaces(t_def, CT)

            # evaluate full constraint Table with less surfaces involved
            CT = build_c_table_from_solids(
                Box=box,
                SurfInfo=(tuple(new_def.get_surfaces_numbers()), Surfs),
                options=options,
                option="full",
            )

            if new_def.operator == "AND":
                new_def.simplify(CT)
                new_def.clean()
            else:
                for i, s in enumerate(simplify):
                    if not s:
                        continue
                    comp = new_def.elements[i]
                    comp.simplify(CT)
                    comp.clean()

            m.set_definition(new_def)
            m.Definition.join_operators()
            m.Definition.level_update()


# TODO this function looks like it is not used in the code.
def extra_plane_cyl_face(face, box, surfaces):
    wire = face.OuterWire
    planes_id = []
    for e in wire.OrderedEdges:
        curve = str(e.Curve)
        if curve[0:6] == "Circle" or curve == "<Ellipse object>":
            dir = e.Curve.Axis
            center = e.Curve.Center
            if curve == "<Ellipse object>":
                dim1 = e.Curve.MinorRadius
                dim2 = e.Curve.MajorRadius
            else:
                dim1 = e.Curve.Radius
                dim2 = e.Curve.Radius
            plane = GeounedSurface(
                ("Plane", (center, dir, dim1, dim2)), box, Face="Build"
            )
            id, exist = surfaces.add_plane(plane)
            if exist:
                pp = surfaces.get_surface(id)
                if is_opposite(plane.Surf.Axis, pp.Surf.Axis, pln_angle):
                    id = -id
            planes_id.append(id)
    return planes_id


def add_cone_plane(definition, cones_list, surfaces, universe_box, angle):
    x_axis = FreeCAD.Vector(1, 0, 0)
    y_axis = FreeCAD.Vector(0, 1, 0)
    z_axis = FreeCAD.Vector(0, 0, 1)

    for cid in cones_list:
        cone = surfaces.get_surface(abs(cid))
        if (
            is_parallel(vector_1=cone.Surf.Axis, vector_2=x_axis, tolerance=angle)
            or is_parallel(vector_1=cone.Surf.Axis, vector_2=y_axis, tolerance=angle)
            or is_parallel(vector_1=cone.Surf.Axis, vector_2=z_axis, tolerance=angle)
        ):
            continue

        plane = GeounedSurface(
            ("Plane", (cone.Surf.Apex, cone.Surf.Axis, 1, 1)),
            universe_box,
            Face="Build",
        )
        pid, exist = surfaces.add_plane(plane)

        if exist:
            p = surfaces.get_surface(pid)
            if is_opposite(
                vector_1=plane.Surf.Axis, vector_2=p.Surf.Axis, tolerance=pln_angle
            ):
                pid = -pid

        if cid > 0:
            insert_in_sequence(definition, cid, -pid, "OR")
        else:
            insert_in_sequence(definition, cid, pid, "AND")
