import logging

import FreeCAD
import Part

from ..decompose import decom_one as Decom
from ..utils import basic_functions_part2 as BF
from ..utils import geometry_gu as GU
from ..utils.basic_functions_part1 import is_opposite, is_parallel
from ..utils.boolean_function import BoolSequence

logger = logging.getLogger("general_logger")


def commonEdge(face1, face2):
    for e1 in face1.Edges:
        for e2 in face2.Edges:
            if e1.isSame(e2):
                return e1
    return None


def isConvex(face1, face2, edge):
    de = 0.1
    tol = 1.0e-5
    delta_point = edge.Vertexes[1].Point - edge.Vertexes[0].Point
    delta_point.normalize()
    axis = delta_point.cross(face2.Surface.Axis)

    point = edge.CenterOfMass + de * axis
    if not face2.__face__.isInside(point, tol, True):
        axis = -axis

    convex = False
    if face1.Surface.Axis.dot(axis) < 0:
        if face1.Orientation == "Forward":
            convex = True
    else:
        if face1.Orientation == "Reversed":
            convex = True
    return convex


def removeElement(Faces, idf):
    for i, f in enumerate(Faces):
        if f[0] == idf:
            del Faces[i]
            break


def is_inverted(solid):
    face = solid.Faces[0]
    parameter_range = face.ParameterRange
    u = (parameter_range[1] + parameter_range[0]) / 2.0
    v = (parameter_range[3] + parameter_range[2]) / 2.0

    point2 = face.CenterOfMass.add(face.normalAt(u, v).multiply(1.0e-6))

    if solid.isInside(point2, 1e-7, False):
        return True
    else:
        return False


# TODO rename this function as there are two with the same name
def get_id(face_in, Surfaces, options, tolerances, numeric_format):

    if is_parallel(face_in.Axis, FreeCAD.Vector(1, 0, 0), tolerances.pln_angle):
        plane = "PX"
    elif is_parallel(face_in.Axis, FreeCAD.Vector(0, 1, 0), tolerances.pln_angle):
        plane = "PY"
    elif is_parallel(face_in.Axis, FreeCAD.Vector(0, 0, 1), tolerances.pln_angle):
        plane = "PZ"
    else:
        plane = "P"

    for s in Surfaces[plane]:
        if BF.is_same_plane(
            face_in,
            s.Surf,
            options=options,
            tolerances=tolerances,
            numeric_format=numeric_format,
        ):
            return s.Index

    return 0


def translate(meta_list, surfaces, universe_box, setting, options, tolerances, numeric_format):
    tot_solid = len(meta_list)
    for i, m in enumerate(meta_list):
        if m.IsEnclosure:
            continue
        logger.info(f"Decomposing solid: {i}/{tot_solid}")
        if setting.debug:
            logger.info(m.Comments)
            if m.IsEnclosure:
                m.Solids[0].exportStep(f"origEnclosure_{i}.stp")
            else:
                m.Solids[0].exportStep(f"origSolid_{i}.stp")

        surfaces.extend(
            Decom.extract_surfaces(
                Part.makeCompound(m.Solids),
                "Plane3Pts",
                universe_box,
                options,
                tolerances,
                numeric_format,
                MakeObj=False,
            )
        )
        set_definition(m, surfaces, options, tolerances)


# TODO rename this, but be careful as there are other functions in the code with the same name
def set_definition(meta_obj, surfaces, options, tolerances):
    solids = meta_obj.Solids
    s_def = BoolSequence(operator="OR")

    for sol in solids:
        subSol = BoolSequence(operator="AND")
        flag_inv = is_inverted(sol)
        solid_gu = GU.SolidGu(sol, tolerances=tolerances, plane3Pts=True)

        faces = []
        for face in solid_gu.Faces:
            if abs(face.Area) < 1e-2:
                continue
            if face.Area < 0:
                logger.warning("Negative surface Area")
            if str(face.Surface) != "<Plane object>":
                logger.warning("All surfaces must be plane")
                continue
            if face.Orientation not in ("Forward", "Reversed"):
                continue

            id = get_id(face.Surface, surfaces, options, tolerances)
            s = surfaces.get_surface(id)
            if is_opposite(face.Surface.Axis, s.Surf.Axis, tolerances.pln_angle):
                id = -id
            if face.Orientation == "Forward":
                id = -id
            if flag_inv:
                id = -id
            faces.append((id, face))

        while len(faces) > 0:
            id1, face1 = faces[0]
            no_convex = []
            for id2, face2 in faces[1:]:
                edge = commonEdge(face1, face2)
                if edge is None:
                    continue
                if not isConvex(face1, face2, edge):
                    no_convex.append(id2)

            if no_convex != []:
                no_convex.insert(0, id1)
                for i in no_convex:
                    removeElement(faces, i)
                orPlanes = BoolSequence(operator="OR")
                orPlanes.append(*no_convex)
                subSol.append(orPlanes)
            else:
                removeElement(faces, id1)
                subSol.append(id1)

        s_def.append(subSol)

    meta_obj.set_definition(s_def)
