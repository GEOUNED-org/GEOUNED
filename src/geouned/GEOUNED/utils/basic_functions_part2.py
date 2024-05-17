#
# Set of useful functions used in different parts of the code
#
import logging
import math

import FreeCAD

from .basic_functions_part1 import (
    is_in_line,
    is_in_plane,
    is_in_tolerance,
    is_opposite,
    is_parallel,
    is_same_value,
)
from ..write.functions import mcnp_surface


def Fuzzy(index, dtype, surf1, surf2, val, tol, options, tolerances, numeric_format):

    fuzzy_logger = logging.getLogger("fuzzy_logger")

    same = val <= tol

    if dtype == "plane":
        p1str = mcnp_surface(index, "Plane", surf1, options, tolerances, numeric_format)
        p2str = mcnp_surface(0, "Plane", surf2, options, tolerances, numeric_format)
        fuzzy_logger.info(f"Same surface : {same}")
        fuzzy_logger.info(f"Plane distance / Tolerance : {val} {tol}\n {p1str}\n {p2str}\n\n")

    elif dtype == "cylRad":
        cyl1str = mcnp_surface(index, "Cylinder", surf1, options, tolerances, numeric_format)
        cyl2str = mcnp_surface(0, "Cylinder", surf2, options, tolerances, numeric_format)
        fuzzy_logger.info(f"Same surface : {same}")
        fuzzy_logger.info(f"Diff Radius / Tolerance: {val} {tol}")
        fuzzy_logger.info(f"{cyl1str}\n {cyl2str}\n\n")

    elif dtype == "cylAxs":
        cyl1str = mcnp_surface(index, "Cylinder", surf1, options, tolerances, numeric_format)
        cyl2str = mcnp_surface(0, "Cylinder", surf2, options, tolerances, numeric_format)
        fuzzy_logger.info(f"Same surface : {same}")
        fuzzy_logger.info(f"Dist Axis / Tolerance: {val} {tol}")
        fuzzy_logger.info(f"{cyl1str} {cyl2str}")


def is_same_plane(p1, p2, options, tolerances, numeric_format, fuzzy=(False, 0)):
    if is_parallel(p1.Axis, p2.Axis, tolerances.pln_angle):
        d1 = p1.Axis.dot(p1.Position)
        d2 = p2.Axis.dot(p2.Position)
        if is_opposite(p1.Axis, p2.Axis, tolerances.pln_angle):
            d2 = -d2
        d = abs(d1 - d2)
        if tolerances.relativeTol:
            tol = tolerances.pln_distance * max(p2.dimL1, p2.dimL2)
        else:
            tol = tolerances.pln_distance

        isSame, is_fuzzy = is_in_tolerance(d, tol, 0.5 * tol, 2 * tol)
        if is_fuzzy and fuzzy[0]:
            Fuzzy(fuzzy[1], "plane", p2, p1, d, tol, options, tolerances, numeric_format)
        return isSame
    return False


def is_same_cylinder(
    cyl1,
    cyl2,
    options,
    tolerances,
    numeric_format,
    fuzzy=(False, 0),
):
    if tolerances.relativeTol:
        rtol = tolerances.cyl_distance * max(cyl2.Radius, cyl1.Radius)
    else:
        rtol = tolerances.cyl_distance

    is_same_rad, is_fuzzy = is_in_tolerance(cyl2.Radius - cyl1.Radius, rtol, 0.5 * rtol, 2 * rtol)
    if is_fuzzy and fuzzy[0]:
        Fuzzy(
            fuzzy[1],
            "cylRad",
            cyl2,
            cyl1,
            abs(cyl2.Radius - cyl1.Radius),
            rtol,
            options,
            tolerances,
            numeric_format,
        )

    if is_same_rad:
        if is_parallel(cyl1.Axis, cyl2.Axis, tolerances.cyl_angle):
            c12 = cyl1.Center - cyl2.Center
            d = cyl1.Axis.cross(c12).Length

            if tolerances.relativeTol:
                tol = tolerances.cyl_distance * max(cyl1.Center.Length, cyl2.Center.Length)
            else:
                tol = tolerances.cyl_distance

            is_same_center, is_fuzzy = is_in_tolerance(d, tol, 0.5 * tol, 2 * tol)
            if is_fuzzy and fuzzy[0]:
                Fuzzy(
                    fuzzy[1],
                    "cylAxs",
                    cyl1,
                    cyl2,
                    d,
                    tol,
                    options,
                    tolerances,
                    numeric_format,
                )

            return is_same_center
    return False


def is_same_cone(cone1, cone2, dtol=1e-6, atol=1e-6, rel_tol=True):
    if is_same_value(cone1.SemiAngle, cone2.SemiAngle, atol):
        if is_parallel(cone1.Axis, cone2.Axis, atol):
            if rel_tol:
                tol = dtol * max(cone1.Apex.Length, cone2.Apex.Length)
            else:
                tol = dtol
            return cone1.Apex.isEqual(cone2.Apex, tol)
    return False


def is_same_sphere(sph1, sph2, tolerance=1e-6, rel_tol=True):
    if rel_tol:
        rtol = tolerance * max(sph2.Radius, sph1.Radius)
    else:
        rtol = tolerance
    if is_same_value(sph1.Radius, sph2.Radius, rtol):
        if rel_tol:
            ctol = tolerance * max(sph2.Center.Length, sph1.Center.Length)
        else:
            ctol = tolerance
        return sph1.Center.isEqual(sph2.Center, ctol)

    return False


def is_same_torus(tor1, tor2, dtol=1e-6, atol=1e-6, rel_tol=True):
    if is_parallel(tor1.Axis, tor2.Axis, atol):
        if tor1.Axis.dot(tor2.Axis) < 0:
            return False  # Assume same cone with oposite axis as different
        if rel_tol:
            Rtol = dtol * max(tor1.MajorRadius, tor2.MajorRadius)
            rtol = dtol * max(tor1.MinorRadius, tor2.MinorRadius)
        else:
            Rtol = dtol
            rtol = dtol

        if is_same_value(tor1.MajorRadius, tor2.MajorRadius, Rtol) and is_same_value(tor1.MinorRadius, tor2.MinorRadius, rtol):
            if rel_tol:
                ctol = dtol * max(tor1.Center.Length, tor2.Center.Length)
            else:
                ctol = dtol
            return tor1.Center.isEqual(tor2.Center, ctol)
    return False


def is_duplicate_in_list(num_str1, i, lista):
    for j, elem2 in enumerate(lista):
        if i == j:
            continue
        num_str2 = f"{elem2:11.4E}"
        num_str3 = f"{elem2 + 2.0 * math.pi:11.4E}"
        num_str4 = f"{elem2 - 2.0 * math.pi:11.4E}"

        if abs(float(num_str2)) < 1.0e-5:
            num_str2 = "%11.4E" % 0.0

        if abs(float(num_str3)) < 1.0e-5:
            num_str3 = "%11.4E" % 0.0

        if abs(float(num_str4)) < 1.0e-5:
            num_str4 = "%11.4E" % 0.0

        if num_str1 == num_str2 or num_str1 == num_str3 or num_str1 == num_str4:
            return True

    return False
