import FreeCAD
import Part

from ..Decompose import Decom_one as Decom
from ..Utils import BasicFunctions_part2 as BF
from ..Utils import Geometry_GU as GU
from ..Utils.BasicFunctions_part1 import isOposite, isParallel
from ..Utils.booleanFunction import BoolSequence
from ..Utils.Options.Classes import Options as opt
from ..Utils.Options.Classes import Tolerances as tol


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


def isInverted(solid):
    face = solid.Faces[0]
    parameter_range = face.ParameterRange
    u = (parameter_range[1] + parameter_range[0]) / 2.0
    v = (parameter_range[3] + parameter_range[2]) / 2.0

    point2 = face.CenterOfMass.add(face.normalAt(u, v).multiply(1.0e-6))

    if solid.isInside(point2, 1e-7, False):
        return True
    else:
        return False


def getId(face_in, Surfaces):

    if isParallel(face_in.Axis, FreeCAD.Vector(1, 0, 0), tol.pln_angle):
        plane = "PX"
    elif isParallel(face_in.Axis, FreeCAD.Vector(0, 1, 0), tol.pln_angle):
        plane = "PY"
    elif isParallel(face_in.Axis, FreeCAD.Vector(0, 0, 1), tol.pln_angle):
        plane = "PZ"
    else:
        plane = "P"

    for s in Surfaces[plane]:
        if BF.isSamePlane(
            face_in,
            s.Surf,
            dtol=tol.pln_distance,
            atol=tol.pln_angle,
            rel_tol=tol.relativeTol,
        ):
            return s.Index

    return 0


def translate(meta_list, surfaces, universe_box, setting):
    tot_solid = len(meta_list)
    for i, m in enumerate(meta_list):
        if m.IsEnclosure:
            continue
        print(f"Decomposing solid: {i}/{tot_solid} ")
        if setting["debug"]:
            print(m.Comments)
            if m.IsEnclosure:
                m.Solids[0].exportStep(f"origEnclosure_{i}.stp")
            else:
                m.Solids[0].exportStep(f"origSolid_{i}.stp")

        surfaces.extend(
            Decom.ExtractSurfaces(
                Part.makeCompound(m.Solids), "Plane3Pts", universe_box, MakeObj=False
            )
        )
        setDefinition(m, surfaces)


def setDefinition(meta_obj, surfaces):
    solids = meta_obj.Solids
    s_def = BoolSequence(operator="OR")

    for sol in solids:
        subSol = BoolSequence(operator="AND")
        flag_inv = isInverted(sol)
        solid_gu = GU.SolidGu(sol, plane3Pts=True)

        faces = []
        for face in solid_gu.Faces:
            if abs(face.Area) < 1e-2:
                continue
            if face.Area < 0:
                if opt.verbose:
                    print("Warning : Negative surface Area")
            if str(face.Surface) != "<Plane object>":
                print("Warning : All surfaces must be plane")
                continue
            if face.Orientation not in ("Forward", "Reversed"):
                continue

            id = getId(face.Surface, surfaces)
            s = surfaces.getSurface(id)
            if isOposite(face.Surface.Axis, s.Surf.Axis, tol.pln_angle):
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

    meta_obj.setDefinition(s_def)
