#
# Set of useful functions used in different parts of the code
#
import math

import FreeCAD

from ..Utils.BasicFunctions_part1 import (
    isInLine,
    isInPlane,
    isInTolerance,
    isOposite,
    isParallel,
    isSameValue,
)
from ..Write.Functions import MCNPSurface

same_surf_fic = open("fuzzySurfaces", "w")


def Fuzzy(index, dtype, surf1, surf2, val, tol):

    same = val <= tol

    if dtype == "plane":
        p1str = MCNPSurface(index, "Plane", surf1)
        p2str = MCNPSurface(0, "Plane", surf2)
        line = "Same surface : {}\nPlane distance / Tolerance : {} {}\n {}\n {}\n\n".format(
            same, val, tol, p1str, p2str
        )
        same_surf_fic.write(line)

    elif dtype == "cylRad":
        cyl1str = MCNPSurface(index, "Cylinder", surf1)
        cyl2str = MCNPSurface(0, "Cylinder", surf2)
        line = "Same surface : {}\nDiff Radius / Tolerance: {} {}\n {}\n {}\n\n".format(
            same, val, tol, cyl1str, cyl2str
        )
        same_surf_fic.write(line)

    elif dtype == "cylAxs":
        cyl1str = MCNPSurface(index, "Cylinder", surf1)
        cyl2str = MCNPSurface(0, "Cylinder", surf2)
        line = "Same surface : {}\nDist Axis / Tolerance: {} {}\n {}\n {}\n\n".format(
            same, val, tol, cyl1str, cyl2str
        )
        same_surf_fic.write(line)


def isSamePlane(p1, p2, dtol=1e-6, atol=1e-6, rel_tol=True, fuzzy=(False, 0)):
    if isParallel(p1.Axis, p2.Axis, atol):
        d1 = p1.Axis.dot(p1.Position)
        d2 = p2.Axis.dot(p2.Position)
        if isOposite(p1.Axis, p2.Axis, atol):
            d2 = -d2
        d = abs(d1 - d2)
        if rel_tol:
            tol = dtol * max(p2.dim1, p2.dim2)
        else:
            tol = dtol

        isSame, is_fuzzy = isInTolerance(d, tol, 0.5 * tol, 2 * tol)
        if is_fuzzy and fuzzy[0]:
            Fuzzy(fuzzy[1], "plane", p2, p1, d, tol)
        return isSame
    return False


def isSameCylinder(cyl1, cyl2, dtol=1e-6, atol=1e-6, rel_tol=True, fuzzy=(False, 0)):
    if rel_tol:
        rtol = dtol * max(cyl2.Radius, cyl1.Radius)
    else:
        rtol = dtol

    is_same_rad, is_fuzzy = isInTolerance(
        cyl2.Radius - cyl1.Radius, rtol, 0.5 * rtol, 2 * rtol
    )
    if is_fuzzy and fuzzy[0]:
        Fuzzy(fuzzy[1], "cylRad", cyl2, cyl1, abs(cyl2.Radius - cyl1.Radius), rtol)

    if is_same_rad:
        if isParallel(cyl1.Axis, cyl2.Axis, atol):
            c12 = cyl1.Center - cyl2.Center
            d = cyl1.Axis.cross(c12).Length

            if rel_tol:
                tol = dtol * max(cyl1.Center.Length, cyl2.Center.Length)
            else:
                tol = dtol

            is_same_center, is_fuzzy = isInTolerance(d, tol, 0.5 * tol, 2 * tol)
            if is_fuzzy and fuzzy[0]:
                Fuzzy(fuzzy[1], "cylAxs", cyl1, cyl2, d, tol)

            return is_same_center
    return False


def isSameCone(cone1, cone2, dtol=1e-6, atol=1e-6, rel_tol=True):
    if isSameValue(cone1.SemiAngle, cone2.SemiAngle, atol):
        if isParallel(cone1.Axis, cone2.Axis, atol):
            if rel_tol:
                tol = dtol * max(cone1.Apex.Length, cone2.Apex.Length)
            else:
                tol = dtol
            return cone1.Apex.isEqual(cone2.Apex, tol)
    return False


def isSameSphere(sph1, sph2, tolerance=1e-6, rel_tol=True):
    if rel_tol:
        rtol = tolerance * max(sph2.Radius, sph1.Radius)
    else:
        rtol = tolerance
    if isSameValue(sph1.Radius, sph2.Radius, rtol):
        if rel_tol:
            ctol = tolerance * max(sph2.Center.Length, sph1.Center.Length)
        else:
            ctol = tolerance
        return sph1.Center.isEqual(sph2.Center, ctol)

    return False


def isSameTorus(tor1, tor2, dtol=1e-6, atol=1e-6, rel_tol=True):
    if isParallel(tor1.Axis, tor2.Axis, atol):
        if tor1.Axis.dot(tor2.Axis) < 0:
            return False  # Assume same cone with oposite axis as different
        if rel_tol:
            Rtol = dtol * max(tor1.MajorRadius, tor2.MajorRadius)
            rtol = dtol * max(tor1.MinorRadius, tor2.MinorRadius)
        else:
            Rtol = dtol
            rtol = dtol

        if isSameValue(tor1.MajorRadius, tor2.MajorRadius, Rtol) and isSameValue(
            tor1.MinorRadius, tor2.MinorRadius, rtol
        ):
            if rel_tol:
                ctol = dtol * max(tor1.Center.Length, tor2.Center.Length)
            else:
                ctol = dtol
            return tor1.Center.isEqual(tor2.Center, ctol)
    return False


def isDuplicateInList(num_str1, i, lista):
    for j, elem2 in enumerate(lista):
        if i == j:
            continue
        num_str2 = "%11.4E" % (elem2)
        num_str3 = "%11.4E" % (elem2 + 2.0 * math.pi)
        num_str4 = "%11.4E" % (elem2 - 2.0 * math.pi)

        if abs(float(num_str2)) < 1.0e-5:
            num_str2 = "%11.4E" % 0.0

        if abs(float(num_str3)) < 1.0e-5:
            num_str3 = "%11.4E" % 0.0

        if abs(float(num_str4)) < 1.0e-5:
            num_str4 = "%11.4E" % 0.0

        if num_str1 == num_str2 or num_str1 == num_str3 or num_str1 == num_str4:
            return True

    return False


def isInFaces(face, faces):

    if faces == []:
        return False
    vector_nulo = FreeCAD.Vector(0, 0, 0)
    surface = face.Surface
    kind_surf = str(face.Surface)
    if kind_surf == "<Plane object>":
        axis = surface.Axis
        position = surface.Position

    elif kind_surf == "<Cylinder object>":
        axis = surface.Axis
        radius = surface.Radius
        center = surface.Center

    elif kind_surf == "<Cone object>":
        axis = surface.Axis
        apex = surface.Apex
        semi_angle = surface.SemiAngle

    elif kind_surf[0:6] == "Sphere":
        center = surface.Center
        radius = surface.Radius

    elif kind_surf == "<Toroid object>":
        axis = surface.Axis
        center = surface.Center
        major_radius = surface.MajorRadius
        minor_radius = surface.MinorRadius

    for elem in faces:
        surf = elem.Surface
        ##print surf
        if str(surf) == "<Plane object>" and kind_surf == "<Plane object>":
            vector_cross = axis.cross(surf.Axis)
            if vector_cross == vector_nulo and isInPlane(position, surf):
                return True
        elif str(surf) == "<Cylinder object>" and kind_surf == "<Cylinder object>":
            dir = surf.Axis
            rad = surf.Radius
            pnt = surf.Center
            vector_cross = axis.cross(surf.Axis)
            if (
                vector_cross == vector_nulo
                and radius == rad
                and isInLine(center, dir, pnt)
            ):
                return True
        elif str(surf) == "<Cone object>" and kind_surf == "<Cone object>":
            # corresponding logic for cone
            dir = surf.Axis
            punta = surf.Apex
            semiangle = surf.SemiAngle
            if (
                axis.isEqual(dir, 1e-5)
                and apex.isEqual(punta, 1e-5)
                and (semi_angle - semiangle) < 1e-6
            ):
                return True
        elif str(surf)[0:6] == "Sphere" and kind_surf[0:6] == "Sphere":
            # corresponding logic for sphere
            rad = surf.Radius
            pnt = surf.Center
            if center == pnt and radius == rad:
                return True

        elif str(surf) == "<Toroid object>" and kind_surf == "<Toroid object>":
            # corresponding logic for Torus
            rad_maj = surf.MajorRadius
            rad_min = surf.MinorRadius
            pnt = surf.Center
            dir = surf.Axis
            if (
                axis.isEqual(dir, 1e-5)
                and center.isEqual(pnt, 1e-5)
                and (major_radius - rad_maj) < 1e-6
            ) and (minor_radius - rad_min) < 1e-6:
                return True

    return False


def isInFaces2(face, faces):

    if faces == []:
        return False
    vector_nulo = FreeCAD.Vector(0, 0, 0)
    surface = face
    kind_surf = face.type
    if kind_surf == "<Plane object>":
        axis = surface.Axis
        position = surface.Position

    elif kind_surf == "<Cylinder object>":
        axis = surface.Axis
        radius = surface.Radius
        center = surface.Center

    elif kind_surf == "<Cone object>":
        axis = surface.Axis
        apex = surface.Apex
        semi_angle = surface.SemiAngle

    elif kind_surf[0:6] == "Sphere":
        center = surface.Center
        radius = surface.Radius

    elif kind_surf == "<Toroid object>":
        axis = surface.Axis
        center = surface.Center
        major_radius = surface.MajorRadius
        minor_radius = surface.MinorRadius

    for elem in faces:
        ##print surf
        if elem.type == "<Plane object>" and kind_surf == "<Plane object>":
            vector_cross = axis.cross(elem.Axis)
            if vector_cross == vector_nulo and isInPlane(position, elem.Surface):
                # if (isParallel(elem.Axis,elem.Surface.Axis) and isInPlane(Position,elem.Surface)):
                return True
        elif elem.type == "<Cylinder object>" and kind_surf == "<Cylinder object>":
            dir = elem.Axis
            rad = elem.Radius
            pnt = elem.Center
            vector_cross = axis.cross(elem.Axis)
            if (
                vector_cross == vector_nulo
                and radius == rad
                and isInLine(center, dir, pnt)
            ):
                return True
        elif elem.type == "<Cone object>" and kind_surf == "<Cone object>":
            # corresponding logic for cone
            dir = elem.Axis
            punta = elem.Apex
            semiangle = elem.SemiAngle
            if (
                axis.isEqual(dir, 1e-5)
                and apex.isEqual(punta, 1e-5)
                and (semi_angle - semiangle) < 1e-6
            ):
                return True
        elif elem.type == "Sphere" and kind_surf == "Sphere":
            # corresponding logic for sphere
            rad = elem.Radius
            pnt = elem.Center
            if center == pnt and radius == rad:
                return True
        elif elem.type == "<Toroid object>" and kind_surf == "<Toroid object>":
            # corresponding logic for Torus
            rad_maj = elem.MajorRadius
            rad_min = elem.MinorRadius
            pnt = elem.Center
            dir = elem.Axis
            if (
                axis.isEqual(dir, 1e-5)
                and center.isEqual(pnt, 1e-5)
                and (major_radius - rad_maj) < 1e-6
            ) and (minor_radius - rad_min) < 1e-6:
                return True
    return False
