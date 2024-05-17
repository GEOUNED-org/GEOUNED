#
# Script to obtain the Qform of a Cylinder
#
#
#

import math

import FreeCAD


def rotation_matrix(u, v):
    """Definition of the rotation matrix for two vectors"""

    # defintion of the axis of rotation

    Axis = u.cross(v).normalize()

    u.normalize()
    v.normalize()

    cose = u.dot(v)
    seno = u.cross(v).Length

    R = FreeCAD.Matrix()

    # cose=math.cos(angle)
    # seno=math.sin(angle)
    onecos = 1.0 - cose

    # 1st row
    R.A11 = cose + Axis.x**2 * onecos
    R.A12 = Axis.x * Axis.y * onecos - Axis.z * seno
    R.A13 = Axis.x * Axis.z * onecos + Axis.y * seno

    # 2nd row
    R.A21 = Axis.x * Axis.y * onecos + Axis.z * seno
    R.A22 = cose + Axis.y**2 * onecos
    R.A23 = Axis.y * Axis.z * onecos - Axis.x * seno

    # 3rd row
    R.A31 = Axis.z * Axis.x * onecos - Axis.y * seno
    R.A32 = Axis.z * Axis.y * onecos + Axis.x * seno
    R.A33 = cose + Axis.z**2 * onecos

    return R


def q_form_cyl(Axis, Pos, rad):

    R = rotation_matrix(FreeCAD.Vector(1, 0, 0), Axis)
    R.transpose()
    Pos2 = R.multiply(Pos).negative()

    A = R.A21**2 + R.A31**2
    B = R.A22**2 + R.A32**2
    C = R.A23**2 + R.A33**2

    D = 2.0 * (R.A21 * R.A22 + R.A31 * R.A32)
    E = 2.0 * (R.A22 * R.A23 + R.A32 * R.A33)
    F = 2.0 * (R.A23 * R.A21 + R.A33 * R.A31)

    G = 2.0 * (Pos2.y * R.A21 + Pos2.z * R.A31)
    H = 2.0 * (Pos2.y * R.A22 + Pos2.z * R.A32)
    J = 2.0 * (Pos2.y * R.A23 + Pos2.z * R.A33)

    K = Pos2.y**2 + Pos2.z**2 - rad**2

    return (A, B, C, D, E, F, G, H, J, K)


def q_form_cone(Axis, Pos, tan):

    R = rotation_matrix(FreeCAD.Vector(1, 0, 0), Axis)
    R.transpose()
    Pos2 = R.multiply(Pos).negative()

    A = R.A21**2 + R.A31**2 - (tan * R.A11) ** 2
    B = R.A22**2 + R.A32**2 - (tan * R.A12) ** 2
    C = R.A23**2 + R.A33**2 - (tan * R.A13) ** 2

    D = 2.0 * (R.A21 * R.A22 + R.A31 * R.A32 - tan**2 * R.A11 * R.A12)
    E = 2.0 * (R.A22 * R.A23 + R.A32 * R.A33 - tan**2 * R.A12 * R.A13)
    F = 2.0 * (R.A23 * R.A21 + R.A33 * R.A31 - tan**2 * R.A13 * R.A11)

    G = 2.0 * (Pos2.y * R.A21 + Pos2.z * R.A31 - tan**2 * Pos2.x * R.A11)
    H = 2.0 * (Pos2.y * R.A22 + Pos2.z * R.A32 - tan**2 * Pos2.x * R.A12)
    J = 2.0 * (Pos2.y * R.A23 + Pos2.z * R.A33 - tan**2 * Pos2.x * R.A13)

    K = Pos2.y**2 + Pos2.z**2 - (tan * Pos2.x) ** 2

    return (A, B, C, D, E, F, G, H, J, K)
