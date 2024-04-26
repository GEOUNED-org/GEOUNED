import math
import re

import FreeCAD

from ..Utils import Qform as Qform
from ..Utils.BasicFunctions_part1 import isOposite, isParallel
from ..Utils.Options.Classes import McnpNumericFormat as nf
from ..Utils.Options.Classes import Options as opt
from ..Utils.Options.Classes import Tolerances as tol
from .StringFunctions import remove_redundant


class CardLine:
    def __init__(self, card, linesize=80, tabspace=6, fmt=""):
        self.str = card + " "
        self.lineSize = linesize
        self.__leftspace__ = linesize - len(self.str)
        self.tabspace = tabspace
        self.tabstr = " " * tabspace
        self.fmt = fmt

    def extend(self, dataList):
        line = ""
        for data in dataList:
            data_str = "{:{fmt}} ".format(data, fmt=self.fmt)
            if self.__leftspace__ - len(data_str) > 0:
                line += data_str
                self.__leftspace__ -= len(data_str)
            else:
                self.str += line + "\n"
                line = "{}{}".format(self.tabstr, data_str)
                self.__leftspace__ = self.lineSize - len(line)

        self.str += line


class CellString:
    def __init__(self, linesize=80, tabspace=6):
        self.str = ""
        self.lineSize = linesize
        self.__leftspace__ = linesize
        self.tabspace = tabspace
        self.tabstr = " " * tabspace

    def add(self, string):
        self.str += string

    def delLastChar(self):
        self.str = self.str[0:-1]
        self.__leftspace__ += 1

    def wrapLine(self, offset=0):

        self.str = self.str.strip()
        self.str = re.sub(" +", " ", self.str)
        self.str = re.sub(" *\( *", "(", self.str)
        self.str = re.sub(" *\) *", ")", self.str)
        self.str = re.sub(" *: *", ":", self.str)
        self.str = re.sub("\)\(", ") (", self.str)
        self.str = re.sub("(?P<num>\d)\(", "\g<num> (", self.str)
        self.str = re.sub("\)(?P<num>\d)", ") \g<num>", self.str)
        self.str = re.sub("\)-", ") -", self.str)

        if len(self.str) + offset <= self.lineSize:
            return

        newline = ""
        init = 0
        icount = self.lineSize - 1 - offset
        lenleft = len(self.str)

        while True:
            while self.str[icount] not in (")", ":", " "):
                icount -= 1

            newline += self.str[init : icount + 1]
            lenleft -= icount - init + 1
            init = icount + 1
            icount += self.lineSize - self.tabspace
            if lenleft == 0:
                break
            elif icount < len(self.str):
                newline += "\n{}".format(self.tabstr)
            else:
                newline += "\n{}{}".format(self.tabstr, self.str[init:])
                break

        self.str = newline


def changeSurfSign(surf, Seq):
    if Seq.level == 0:
        for i, e in enumerate(Seq.elements):
            if surf == abs(e):
                Seq.elements[i] = -Seq.elements[i]
    else:
        for i, e in enumerate(Seq.elements):
            if type(e) is int:
                if surf == abs(e):
                    Seq.elements[i] = -Seq.elements[i]
            else:
                changeSurfSign(surf, e)


def writeMCNPCellDef(Definition, tabspace=0, offset=0):
    sdef = CellString(tabspace=tabspace)
    strDef = remove_redundant(writeSequenceMCNP(Definition))
    sdef.add(strDef)
    sdef.wrapLine(offset)
    return sdef.str


def writeSerpentCellDef(Definition, tabspace=0, offset=0):
    sdef = CellString(tabspace=tabspace)
    strDef = remove_redundant(writeSequenceSerpent(Definition))
    sdef.add(strDef)
    sdef.wrapLine(offset)
    return sdef.str


def writePHITSCellDef(Definition, tabspace=0, offset=0):
    sdef = CellString(tabspace=tabspace)
    strDef = remove_redundant(writeSequencePHITS(Definition))
    sdef.add(strDef)
    sdef.wrapLine(offset)
    return sdef.str


def writeOpenMCregion(Definition, Wtype="XML"):
    if Wtype == "XML":
        return writeSequenceOMCXML(Definition)
    if Wtype == "PY":
        return writeSequenceOMCPY(Definition)


def writeSequenceMCNP(Seq):
    if Seq.level == 0:
        if Seq.operator == "AND":
            line = "({})".format(" ".join(map(str, Seq.elements)))
        else:
            line = "({})".format(":".join(map(str, Seq.elements)))
    else:
        terms = []
        for e in Seq.elements:
            if type(e) is int:
                terms.append(str(e))
            else:
                terms.append(writeSequenceMCNP(e))

        if Seq.operator == "AND":
            line = "({})".format(" ".join(terms))
        else:
            line = "({})".format(":".join(terms))

    return line


def writeSequenceSerpent(Seq):
    if Seq.level == 0:
        if Seq.operator == "AND":
            line = "({})".format(" ".join(map(str, Seq.elements)))
        else:
            line = "({})".format(":".join(map(str, Seq.elements)))
    else:
        terms = []
        for e in Seq.elements:
            if type(e) is int:
                terms.append(str(e))
            else:
                terms.append(writeSequenceMCNP(e))

        if Seq.operator == "AND":
            line = "({})".format(" ".join(terms))
        else:
            line = "({})".format(":".join(terms))

    return line


def writeSequencePHITS(Seq):
    if Seq.level == 0:
        if Seq.operator == "AND":
            line = "({})".format(" ".join(map(str, Seq.elements)))
        else:
            line = "({})".format(":".join(map(str, Seq.elements)))
    else:
        terms = []
        for e in Seq.elements:
            if type(e) is int:
                terms.append(str(e))
            else:
                terms.append(writeSequencePHITS(e))

        if Seq.operator == "AND":
            line = "({})".format(" ".join(terms))
        else:
            line = "({})".format(":".join(terms))

    return line


def writeSequenceOMCXML(Seq):
    if Seq.level == 0:
        if Seq.operator == "AND":
            line = "({})".format(" ".join(map(str, Seq.elements)))
        else:
            line = "({})".format(" | ".join(map(str, Seq.elements)))
    else:
        terms = []
        for e in Seq.elements:
            if type(e) is int:
                terms.append(str(e))
            else:
                terms.append(writeSequenceOMCXML(e))

        if Seq.operator == "AND":
            line = "({})".format(" ".join(terms))
        else:
            line = "({})".format(" | ".join(terms))
    return line


def writeSequenceOMCPY(Seq, prefix="S"):

    strSurf = lambda surf: (
        "-{}{}".format(prefix, -surf) if surf < 0 else "+{}{}".format(prefix, surf)
    )

    if Seq.level == 0:
        if Seq.operator == "AND":
            line = "({})".format(" & ".join(map(strSurf, Seq.elements)))
        else:
            line = "({})".format(" | ".join(map(strSurf, Seq.elements)))
    else:
        terms = []
        for e in Seq.elements:
            if type(e) is int:
                terms.append(strSurf(e))
            else:
                terms.append(writeSequenceOMCPY(e))

        if Seq.operator == "AND":
            line = "({})".format(" & ".join(terms))
        else:
            line = "({})".format(" | ".join(terms))
    return line


def MCNPSurface(id, Type, surf):
    MCNP_def = ""

    if Type == "Plane":
        if surf.pointDef and opt.prnt3PPlane:
            P1 = surf.Points[0]
            P2 = surf.Points[1]
            P3 = surf.Points[2]
            MCNP_def = """{:<6d} P   {P1[0]:{d}} {P1[1]:{d}} {P1[2]:{d}} 
           {P2[0]:{d}} {P2[1]:{d}} {P2[2]:{d}}
           {P3[0]:{d}} {P3[1]:{d}} {P3[2]:{d}}""".format(
                id, P1=P1 / 10, P2=P2 / 10, P3=P3 / 10, d=nf.P_d
            )
        else:
            A = surf.Axis.x
            B = surf.Axis.y
            C = surf.Axis.z
            D = surf.Axis.dot(surf.Position)
            if surf.Axis.isEqual(FreeCAD.Vector(1, 0, 0), tol.pln_angle):
                MCNP_def = "{:<6d} PX  {:{x}}".format(id, D / 10.0, x=nf.P_xyz)
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 1, 0), tol.pln_angle):
                MCNP_def = "{:<6d} PY  {:{y}}".format(id, D / 10.0, y=nf.P_xyz)
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 0, 1), tol.pln_angle):
                MCNP_def = "{:<6d} PZ  {:{z}}".format(id, D / 10.0, z=nf.P_xyz)
            else:
                MCNP_def = "{:<6d} P   {:{abc}} {:{abc}} {:{abc}} {:{d}}".format(
                    id, A, B, C, D / 10.0, abc=nf.P_abc, d=nf.P_d
                )

    elif Type == "Cylinder":
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        Pos = surf.Center * 0.1
        rad = surf.Radius * 0.1
        if isParallel(Dir, FreeCAD.Vector(1, 0, 0), tol.angle):
            if Pos.y == 0.0 and Pos.z == 0.0:
                MCNP_def = "{:<6d} CX  {:{r}}".format(id, rad, r=nf.C_r)
            else:
                MCNP_def = "{:<6d} C/X  {:{yz}} {:{yz}} {:{r}}".format(
                    id, Pos.y, Pos.z, rad, yz=nf.C_xyz, r=nf.C_r
                )
        elif isParallel(Dir, FreeCAD.Vector(0, 1, 0), tol.angle):
            if Pos.x == 0.0 and Pos.z == 0.0:
                MCNP_def = "{:<6d} CY  {:{r}}".format(id, rad, r=nf.C_r)
            else:
                MCNP_def = "{:<6d} C/Y  {:{xz}} {:{xz}} {:{r}}".format(
                    id, Pos.x, Pos.z, rad, xz=nf.C_xyz, r=nf.C_r
                )
        elif isParallel(Dir, FreeCAD.Vector(0, 0, 1), tol.angle):
            if Pos.y == 0.0 and Pos.x == 0.0:
                MCNP_def = "{:<6d} CZ  {:{r}}".format(id, rad, r=nf.C_r)
            else:
                MCNP_def = "{:<6d} C/Z  {:{xy}} {:{xy}} {:{r}}".format(
                    id, Pos.x, Pos.y, rad, xy=nf.C_xyz, r=nf.C_r
                )
        else:
            # Is not still working fine
            Q = Qform.QFormCyl(Dir, Pos, rad)
            MCNP_def = """\
{:<6d} GQ  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id, v=Q, aTof=nf.GQ_1to6, gToi=nf.GQ_7to9, j=nf.GQ_10
            )

        # Si se quiere rcc en vez de Q form
        #    pnt = surf.Center.sub(surf.Axis.multiply(1.0e7)) # mas alla de 100 m
        #    dir = surf.Axis.multiply(1.0e8)
        #    Vx= pnt.x/10.0
        #    Vy= pnt.y/10.0
        #    Vz= pnt.z/10.0
        #    Hx= dir.x/10.0
        #    Hy= dir.y/10.0
        #    Hz= dir.z/10.0
        #    rad=surf.Radius/10.0
        #    MCNP_def='%i  RCC  %13.7E %13.7E %13.7E %13.7E\n       %13.7E %13.7E %13.7E' %(id,Vx,Vy,Vz,Hx,Hy,Hz,rad)

    elif Type == "Cone":
        Apex = surf.Apex * 0.1
        Dir = surf.Axis * 0.1
        tan = math.tan(surf.SemiAngle)
        X_dir = FreeCAD.Vector(1, 0, 0)
        Y_dir = FreeCAD.Vector(0, 1, 0)
        Z_dir = FreeCAD.Vector(0, 0, 1)
        if isParallel(Dir, X_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, X_dir, tol.angle):
                sheet = -1
            if Apex.y == 0.0 and Apex.z == 0.0:
                MCNP_def = "{:<6d} KX  {:{x}} {:{t2}} {}".format(
                    id, Apex.x, tan**2, sheet, x=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                MCNP_def = "{:<6d} K/X  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=nf.K_xyz,
                    t2=nf.K_tan2,
                )
        elif isParallel(Dir, Y_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, Y_dir, tol.angle):
                sheet = -1
            if Apex.x == 0.0 and Apex.z == 0.0:
                MCNP_def = "{:<6d} KY  {:{y}} {:{t2}} {}".format(
                    id, Apex.y, tan**2, sheet, y=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                MCNP_def = "{:<6d} K/Y  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=nf.K_xyz,
                    t2=nf.K_tan2,
                )
        elif isParallel(Dir, Z_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, Z_dir, tol.angle):
                sheet = -1
            if Apex.x == 0.0 and Apex.y == 0.0:
                MCNP_def = "{:<6d} KZ  {:{z}} {:{t2}} {}".format(
                    id, Apex.z, tan**2, sheet, z=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                MCNP_def = "{:<6d} K/Z  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=nf.K_xyz,
                    t2=nf.K_tan2,
                )
        else:
            Q = Qform.QFormCone(Dir, Apex, tan)
            MCNP_def = """\
{:<6d} GQ  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id, v=Q, aTof=nf.GQ_1to6, gToi=nf.GQ_7to9, j=nf.GQ_10
            )

    elif Type == "Sphere":
        # corresponding logic
        rad = surf.Radius * 0.1
        pnt = surf.Center * 0.1
        if pnt.isEqual(FreeCAD.Vector(0, 0, 0), tol.sph_distance):
            MCNP_def = "{:<6d} SO  {:{r}}".format(id, rad, r=nf.S_r)
        else:
            MCNP_def = "{:<6d} S  {:{xyz}} {:{xyz}} {:{xyz}} {:{r}}".format(
                id, pnt.x, pnt.y, pnt.z, rad, xyz=nf.S_xyz, r=nf.S_r
            )

    elif Type == "Torus":
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        Pos = surf.Center * 0.1
        radMaj = surf.MajorRadius * 0.1
        radMin = surf.MinorRadius * 0.1
        if isParallel(Dir, FreeCAD.Vector(1, 0, 0), tol.angle):
            MCNP_def = """\
{:<6d} TX  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id, Pos.x, Pos.y, Pos.z, radMaj, radMin, radMin, xyz=nf.T_xyz, r=nf.T_r
            )
        elif isParallel(Dir, FreeCAD.Vector(0, 1, 0), tol.angle):
            MCNP_def = """\
{:<6d} TY  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id, Pos.x, Pos.y, Pos.z, radMaj, radMin, radMin, xyz=nf.T_xyz, r=nf.T_r
            )
        elif isParallel(Dir, FreeCAD.Vector(0, 0, 1), tol.angle):
            MCNP_def = """\
{:<6d} TZ  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id, Pos.x, Pos.y, Pos.z, radMaj, radMin, radMin, xyz=nf.T_xyz, r=nf.T_r
            )

    return trim(MCNP_def, 80)


def OpenMCSurface(Type, surf, outXML=True, quadricForm=False):
    if Type == "Plane":
        A = surf.Axis.x
        B = surf.Axis.y
        C = surf.Axis.z
        D = surf.Axis.dot(surf.Position) * 0.1
        if surf.Axis.isEqual(FreeCAD.Vector(1, 0, 0), tol.pln_angle):
            if outXML:
                OMCsurf = "x-plane"
                coeffs = "{:{x}}".format(D, x=nf.P_xyz)
            else:
                OMCsurf = "XPlane"
                coeffs = "x0={}".format(D)

        elif surf.Axis.isEqual(FreeCAD.Vector(0, 1, 0), tol.pln_angle):
            if outXML:
                OMCsurf = "y-plane"
                coeffs = "{:{x}}".format(D, x=nf.P_xyz)
            else:
                OMCsurf = "YPlane"
                coeffs = "y0={}".format(D)

        elif surf.Axis.isEqual(FreeCAD.Vector(0, 0, 1), tol.pln_angle):
            if outXML:
                OMCsurf = "z-plane"
                coeffs = "{:{x}}".format(D, x=nf.P_xyz)
            else:
                OMCsurf = "ZPlane"
                coeffs = "z0={}".format(D)

        else:
            if outXML:
                OMCsurf = "plane"
                coeffs = "{:{abc}} {:{abc}} {:{abc}} {:{d}}".format(
                    A, B, C, D, abc=nf.P_abc, d=nf.P_d
                )
            else:
                OMCsurf = "Plane"
                coeffs = "a={},b={},c={},d={}".format(A, B, C, D)

    elif Type == "Cylinder":
        Pos = surf.Center * 0.1
        Rad = surf.Radius * 0.1
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()

        if isParallel(Dir, FreeCAD.Vector(1, 0, 0), tol.angle):
            if outXML:
                OMCsurf = "x-cylinder"
                coeffs = "{:{xy}} {:{xy}} {:{r}}".format(
                    Pos.y, Pos.z, Rad, xy=nf.C_xyz, r=nf.C_r
                )
            else:
                OMCsurf = "XCylinder"
                coeffs = "y0={},z0={},r={}".format(Pos.y, Pos.z, Rad)

        elif isParallel(Dir, FreeCAD.Vector(0, 1, 0), tol.angle):
            if outXML:
                OMCsurf = "y-cylinder"
                coeffs = "{:{xy}} {:{xy}} {:{r}}".format(
                    Pos.x, Pos.z, Rad, xy=nf.C_xyz, r=nf.C_r
                )
            else:
                OMCsurf = "YCylinder"
                coeffs = "x0={},z0={},r={}".format(Pos.x, Pos.z, Rad)

        elif isParallel(Dir, FreeCAD.Vector(0, 0, 1), tol.angle):
            if outXML:
                OMCsurf = "z-cylinder"
                coeffs = "{:{xy}} {:{xy}} {:{r}}".format(
                    Pos.x, Pos.y, Rad, xy=nf.C_xyz, r=nf.C_r
                )
            else:
                OMCsurf = "ZCylinder"
                coeffs = "x0={},y0={},r={}".format(Pos.x, Pos.y, Rad)

        else:
            if outXML:
                OMCsurf = "quadric"
                Q = Qform.QFormCyl(Dir, Pos, Rad)
                coeffs = "{v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}} {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}} {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}} {v[9]:{j}}".format(
                    v=Q, aTof=nf.GQ_1to6, gToi=nf.GQ_7to9, j=nf.GQ_10
                )
            else:
                if quadricForm:
                    OMCsurf = "Quadric"
                    Q = Qform.QFormCyl(Dir, Pos, Rad)
                    coeffs = "a={v[0]},b={v[1]},c={v[2]},d={v[3]},e={v[4]},f={v[5]},g={v[6]},h={v[7]},j={v[8]},k={v[9]}".format(
                        v=Q
                    )
                else:
                    OMCsurf = "Cylinder"
                    coeffs = "x0={},y0={},z0={},r={},dx={},dy={},dz={}".format(
                        Pos.x, Pos.y, Pos.z, Rad, Dir.x, Dir.y, Dir.z
                    )

    elif Type == "Cone":
        Apex = surf.Apex * 0.1
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        tan = math.tan(surf.SemiAngle)
        tan2 = tan * tan

        X_dir = FreeCAD.Vector(1, 0, 0)
        Y_dir = FreeCAD.Vector(0, 1, 0)
        Z_dir = FreeCAD.Vector(0, 0, 1)

        if isParallel(Dir, X_dir, tol.angle):
            if outXML:
                OMCsurf = "x-cone"
                coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{t2}}".format(
                    Apex.x, Apex.y, Apex.z, tan2, xyz=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                OMCsurf = "XCone"
                coeffs = "x0={},y0={},z0={},r2={}".format(Apex.x, Apex.y, Apex.z, tan2)

        elif isParallel(Dir, Y_dir, tol.angle):
            if outXML:
                OMCsurf = "y-cone"
                coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{t2}}".format(
                    Apex.x, Apex.y, Apex.z, tan2, xyz=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                OMCsurf = "YCone"
                coeffs = "x0={},y0={},z0={},r2={}".format(Apex.x, Apex.y, Apex.z, tan2)

        elif isParallel(Dir, Z_dir, tol.angle):
            if outXML:
                OMCsurf = "z-cone"
                coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{t2}}".format(
                    Apex.x, Apex.y, Apex.z, tan2, xyz=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                OMCsurf = "ZCone"
                coeffs = "x0={},y0={},z0={},r2={}".format(Apex.x, Apex.y, Apex.z, tan2)

        else:
            if outXML:
                OMCsurf = "quadric"
                Q = Qform.QFormCone(Dir, Apex, tan)
                coeffs = "{v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}} {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}} {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}} {v[9]:{j}}".format(
                    v=Q, aTof=nf.GQ_1to6, gToi=nf.GQ_7to9, j=nf.GQ_10
                )
            else:
                if quadricForm:
                    OMCsurf = "Quadric"
                    Q = Qform.QFormCone(Dir, Apex, tan)
                    coeffs = "a={v[0]},b={v[1]},c={v[2]},d={v[3]},e={v[4]},f={v[5]},g={v[6]},h={v[7]},j={v[8]},k={v[9]}".format(
                        v=Q
                    )
                else:
                    OMCsurf = "Cone"
                    coeffs = "x0={},y0={},z0={},r2={},dx={},dy={},dz={}".format(
                        Apex.x, Apex.y, Apex.z, tan2, Dir.x, Dir.y, Dir.z
                    )

    elif Type == "Sphere":
        Center = surf.Center * 0.1
        Rad = surf.Radius * 0.1
        if outXML:
            OMCsurf = "sphere"
            coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{r}}".format(
                Center.x, Center.y, Center.z, Rad, xyz=nf.S_xyz, r=nf.S_r
            )
        else:
            OMCsurf = "Sphere"
            coeffs = "x0={},y0={},z0={},r={}".format(Center.x, Center.y, Center.z, Rad)

    elif Type == "Torus":
        Center = surf.Center * 0.1
        minRad = surf.MinorRadius * 0.1
        majRad = surf.MajorRadius * 0.1
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        if outXML:
            coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{r}} {:{r}} {:{r}}".format(
                Center.x,
                Center.y,
                Center.z,
                majRad,
                minRad,
                minRad,
                xyz=nf.T_xyz,
                r=nf.T_r,
            )
        else:
            coeffs = "x0={},y0={},z0={},r{},r1={},r2={}".format(
                Center.x, Center.y, Center.z, majRad, minRad, minRad
            )

        if isParallel(Dir, FreeCAD.Vector(1, 0, 0), tol.angle):
            OMCsurf = "x-torus" if outXML else "XTorus"
        elif isParallel(Dir, FreeCAD.Vector(0, 1, 0), tol.angle):
            OMCsurf = "y-torus" if outXML else "YTorus"
        elif isParallel(Dir, FreeCAD.Vector(0, 0, 1), tol.angle):
            OMCsurf = "z-torus" if outXML else "ZTorus"
        else:
            OMCsurf = None

    if outXML:
        coeffs = " ".join(coeffs.split())
    return OMCsurf, coeffs


def SerpentSurface(id, Type, surf):
    Serpent_def = ""

    if Type == "Plane":
        if surf.pointDef and opt.prnt3PPlane:
            P1 = surf.Points[0]
            P2 = surf.Points[1]
            P3 = surf.Points[2]
            Serpent_def = f"surf {id} plane {P1.x/10:{nf.P_d}} {P1.y/10:{nf.P_d}} {P1.z/10:{nf.P_d}}\n"
            Serpent_def += (
                f"      {P2.x/10:{nf.P_d}} {P2.y/10:{nf.P_d}} {P2.z/10:{nf.P_d}}\n"
            )
            Serpent_def += (
                f"      {P3.x/10:{nf.P_d}} {P3.y/10:{nf.P_d}} {P3.z/10:{nf.P_d}}"
            )

        else:
            A = surf.Axis.x
            B = surf.Axis.y
            C = surf.Axis.z
            D = surf.Axis.dot(surf.Position)
            if surf.Axis.isEqual(FreeCAD.Vector(1, 0, 0), tol.pln_angle):
                Serpent_def = f"surf {id} px {D/10:{nf.P_xyz}}"
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 1, 0), tol.pln_angle):
                Serpent_def = f"surf {id} py {D/10:{nf.P_xyz}}"
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 0, 1), tol.pln_angle):
                Serpent_def = f"surf {id} pz {D/10:{nf.P_xyz}}"
            else:
                Serpent_def = f"surf {id} plane {A:{nf.P_d}} {B:{nf.P_d}} {C:{nf.P_d}} {D/10:{nf.P_d}}"

    elif Type == "Cylinder":
        Dir = surf.Axis
        Dir.normalize()
        Pos = surf.Center * 0.1
        rad = surf.Radius * 0.1
        if isParallel(Dir, FreeCAD.Vector(1, 0, 0), tol.angle):
            Serpent_def = (
                f"surf {id} cylx {Pos.y:{nf.C_xyz}} {Pos.z:{nf.C_xyz}} {rad:{nf.C_r}}"
            )
        elif isParallel(Dir, FreeCAD.Vector(0, 1, 0), tol.angle):
            Serpent_def = (
                f"surf {id} cyly {Pos.x:{nf.C_xyz}} {Pos.z:{nf.C_xyz}} {rad:{nf.C_r}}"
            )
        elif isParallel(Dir, FreeCAD.Vector(0, 0, 1), tol.angle):
            Serpent_def = f"surf {id} cylz {rad:{nf.C_r}}"
        else:
            # Is not still working fine
            Q = Qform.QFormCyl(Dir, Pos, rad)
            Serpent_def = """\
surf quadratic  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id, v=Q, aTof=nf.GQ_1to6, gToi=nf.GQ_7to9, j=nf.GQ_10
            )

    elif Type == "Cone":
        Apex = surf.Apex * 0.1
        Dir = surf.Axis * 0.1
        tan = math.tan(surf.SemiAngle)
        X_dir = FreeCAD.Vector(1, 0, 0)
        Y_dir = FreeCAD.Vector(0, 1, 0)
        Z_dir = FreeCAD.Vector(0, 0, 1)

        # Need to check this
        # Serpent has no specific card for cone at origin, explicit origin only

        if isParallel(Dir, X_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, X_dir, tol.angle):
                sheet = -1
            Serpent_def = "surf ckx {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                id, Apex.x, Apex.y, Apex.z, tan**2, sheet, xyz=nf.K_xyz, t2=nf.K_tan2
            )
        elif isParallel(Dir, Y_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, Y_dir, tol.angle):
                sheet = -1
            Serpent_def = "surf cky {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                id, Apex.x, Apex.y, Apex.z, tan**2, sheet, xyz=nf.K_xyz, t2=nf.K_tan2
            )
        elif isParallel(Dir, Z_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, Z_dir, tol.angle):
                sheet = -1
            Serpent_def = "surf ckz {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                id, Apex.x, Apex.y, Apex.z, tan**2, sheet, xyz=nf.K_xyz, t2=nf.K_tan2
            )
        else:
            Q = Qform.QFormCone(Dir, Apex, tan)

    elif Type == "Sphere":
        rad = surf.Radius * 0.1
        pnt = surf.Center * 0.1
        # Serpent has only explicit spheres at the origin
        Serpent_def = f"surf {id} sph {pnt.x:{nf.S_xyz}} {pnt.y:{nf.S_xyz}} {pnt.z:{nf.S_xyz}} {rad:{nf.S_r}}"

    elif Type == "Torus":
        Dir = surf.Axis
        Dir.normalize()
        Pos = surf.Center * 0.1
        radMaj = surf.MajorRadius * 0.1
        radMin = surf.MinorRadius * 0.1
        if isParallel(Dir, FreeCAD.Vector(1, 0, 0), tol.angle):
            Serpent_def = f"surf {id} torx {Pos.x:{nf.T_xyz}} {Pos.y:{nf.T_xyz}} {Pos.z:{nf.T_xyz}}\n"
            Serpent_def += (
                f"      {radMaj:{nf.T_r}} {radMin:{nf.T_r}} {radMin:{nf.T_r}}"
            )
        elif isParallel(Dir, FreeCAD.Vector(0, 1, 0), tol.angle):
            Serpent_def = f"surf {id} tory {Pos.x:{nf.T_xyz}} {Pos.y:{nf.T_xyz}} {Pos.z:{nf.T_xyz}}\n"
            Serpent_def += (
                f"      {radMaj:{nf.T_r}} {radMin:{nf.T_r}} {radMin:{nf.T_r}}"
            )
        elif isParallel(Dir, FreeCAD.Vector(0, 0, 1), tol.angle):
            Serpent_def = f"surf {id} torz {Pos.x:{nf.T_xyz}} {Pos.y:{nf.T_xyz}} {Pos.z:{nf.T_xyz}}\n"
            Serpent_def += (
                f"      {radMaj:{nf.T_r}} {radMin:{nf.T_r}} {radMin:{nf.T_r}}"
            )

    return Serpent_def


def PHITSSurface(id, Type, surf):
    PHITS_def = ""

    if Type == "Plane":
        if surf.pointDef and opt.prnt3PPlane:
            P1 = surf.Points[0]
            P2 = surf.Points[1]
            P3 = surf.Points[2]
            PHITS_def = """{:<6d} P   {P1[0]:{d}} {P1[1]:{d}} {P1[2]:{d}} 
           {P2[0]:{d}} {P2[1]:{d}} {P2[2]:{d}}
           {P3[0]:{d}} {P3[1]:{d}} {P3[2]:{d}}""".format(
                id, P1=P1 / 10, P2=P2 / 10, P3=P3 / 10, d=nf.P_d
            )
        else:
            A = surf.Axis.x
            B = surf.Axis.y
            C = surf.Axis.z
            D = surf.Axis.dot(surf.Position)
            if surf.Axis.isEqual(FreeCAD.Vector(1, 0, 0), tol.pln_angle):
                PHITS_def = "{:<6d} PX  {:{x}}".format(id, D / 10.0, x=nf.P_xyz)
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 1, 0), tol.pln_angle):
                PHITS_def = "{:<6d} PY  {:{y}}".format(id, D / 10.0, y=nf.P_xyz)
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 0, 1), tol.pln_angle):
                PHITS_def = "{:<6d} PZ  {:{z}}".format(id, D / 10.0, z=nf.P_xyz)
            else:
                PHITS_def = "{:<6d} P   {:{abc}} {:{abc}} {:{abc}} {:{d}}".format(
                    id, A, B, C, D / 10.0, abc=nf.P_abc, d=nf.P_d
                )

    elif Type == "Cylinder":
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        Pos = surf.Center * 0.1
        rad = surf.Radius * 0.1
        if isParallel(Dir, FreeCAD.Vector(1, 0, 0), tol.angle):
            if Pos.y == 0.0 and Pos.z == 0.0:
                PHITS_def = "{:<6d} CX  {:{r}}".format(id, rad, r=nf.C_r)
            else:
                PHITS_def = "{:<6d} C/X  {:{yz}} {:{yz}} {:{r}}".format(
                    id, Pos.y, Pos.z, rad, yz=nf.C_xyz, r=nf.C_r
                )
        elif isParallel(Dir, FreeCAD.Vector(0, 1, 0), tol.angle):
            if Pos.x == 0.0 and Pos.z == 0.0:
                PHITS_def = "{:<6d} CY  {:{r}}".format(id, rad, r=nf.C_r)
            else:
                PHITS_def = "{:<6d} C/Y  {:{xz}} {:{xz}} {:{r}}".format(
                    id, Pos.x, Pos.z, rad, xz=nf.C_xyz, r=nf.C_r
                )
        elif isParallel(Dir, FreeCAD.Vector(0, 0, 1), tol.angle):
            if Pos.y == 0.0 and Pos.x == 0.0:
                PHITS_def = "{:<6d} CZ  {:{r}}".format(id, rad, r=nf.C_r)
            else:
                PHITS_def = "{:<6d} C/Z  {:{xy}} {:{xy}} {:{r}}".format(
                    id, Pos.x, Pos.y, rad, xy=nf.C_xyz, r=nf.C_r
                )
        else:
            # Is not still working fine
            Q = Qform.QFormCyl(Dir, Pos, rad)
            PHITS_def = """\
{:<6d} GQ  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id, v=Q, aTof=nf.GQ_1to6, gToi=nf.GQ_7to9, j=nf.GQ_10
            )

        # Si se quiere rcc en vez de Q form
        #    pnt = surf.Center.sub(surf.Axis.multiply(1.0e7)) # mas alla de 100 m
        #    dir = surf.Axis.multiply(1.0e8)
        #    Vx= pnt.x/10.0
        #    Vy= pnt.y/10.0
        #    Vz= pnt.z/10.0
        #    Hx= dir.x/10.0
        #    Hy= dir.y/10.0
        #    Hz= dir.z/10.0
        #    rad=surf.Radius/10.0
        #    MCNP_def='%i  RCC  %13.7E %13.7E %13.7E %13.7E\n       %13.7E %13.7E %13.7E' %(id,Vx,Vy,Vz,Hx,Hy,Hz,rad)

    elif Type == "Cone":
        Apex = surf.Apex * 0.1
        Dir = surf.Axis * 0.1
        tan = math.tan(surf.SemiAngle)
        X_dir = FreeCAD.Vector(1, 0, 0)
        Y_dir = FreeCAD.Vector(0, 1, 0)
        Z_dir = FreeCAD.Vector(0, 0, 1)
        if isParallel(Dir, X_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, X_dir, tol.angle):
                sheet = -1
            if Apex.y == 0.0 and Apex.z == 0.0:
                PHITS_def = "{:<6d} KX  {:{x}} {:{t2}} {}".format(
                    id, Apex.x, tan**2, sheet, x=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                PHITS_def = "{:<6d} K/X  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=nf.K_xyz,
                    t2=nf.K_tan2,
                )
        elif isParallel(Dir, Y_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, Y_dir, tol.angle):
                sheet = -1
            if Apex.x == 0.0 and Apex.z == 0.0:
                PHITS_def = "{:<6d} KY  {:{y}} {:{t2}} {}".format(
                    id, Apex.y, tan**2, sheet, y=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                PHITS_def = "{:<6d} K/Y  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=nf.K_xyz,
                    t2=nf.K_tan2,
                )
        elif isParallel(Dir, Z_dir, tol.angle):
            sheet = 1
            if isOposite(Dir, Z_dir, tol.angle):
                sheet = -1
            if Apex.x == 0.0 and Apex.y == 0.0:
                PHITS_def = "{:<6d} KZ  {:{z}} {:{t2}} {}".format(
                    id, Apex.z, tan**2, sheet, z=nf.K_xyz, t2=nf.K_tan2
                )
            else:
                PHITS_def = "{:<6d} K/Z  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=nf.K_xyz,
                    t2=nf.K_tan2,
                )
        else:
            Q = Qform.QFormCone(Dir, Apex, tan)
            PHITS_def = """\
{:<6d} GQ  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id, v=Q, aTof=nf.GQ_1to6, gToi=nf.GQ_7to9, j=nf.GQ_10
            )

    elif Type == "Sphere":
        # corresponding logic
        rad = surf.Radius * 0.1
        pnt = surf.Center * 0.1
        if pnt.isEqual(FreeCAD.Vector(0, 0, 0), tol.sph_distance):
            PHITS_def = "{:<6d} SO  {:{r}}".format(id, rad, r=nf.S_r)
        else:
            PHITS_def = "{:<6d} S  {:{xyz}} {:{xyz}} {:{xyz}} {:{r}}".format(
                id, pnt.x, pnt.y, pnt.z, rad, xyz=nf.S_xyz, r=nf.S_r
            )

    elif Type == "Torus":
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        Pos = surf.Center * 0.1
        radMaj = surf.MajorRadius * 0.1
        radMin = surf.MinorRadius * 0.1
        if isParallel(Dir, FreeCAD.Vector(1, 0, 0), tol.angle):
            PHITS_def = """\
{:<6d} TX  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id, Pos.x, Pos.y, Pos.z, radMaj, radMin, radMin, xyz=nf.T_xyz, r=nf.T_r
            )
        elif isParallel(Dir, FreeCAD.Vector(0, 1, 0), tol.angle):
            PHITS_def = """\
{:<6d} TY  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id, Pos.x, Pos.y, Pos.z, radMaj, radMin, radMin, xyz=nf.T_xyz, r=nf.T_r
            )
        elif isParallel(Dir, FreeCAD.Vector(0, 0, 1), tol.angle):
            PHITS_def = """\
{:<6d} TZ  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id, Pos.x, Pos.y, Pos.z, radMaj, radMin, radMin, xyz=nf.T_xyz, r=nf.T_r
            )

    return trim(PHITS_def, 80)


def trim(surfDef, lineLength=80):

    lines = surfDef.split("\n")
    if len(lines) == 1 and len(lines[0]) <= lineLength:
        return surfDef

    longLine = []
    for i, line in enumerate(lines):
        if len(line) > lineLength:
            longLine.append(i)

    if len(longLine) == 0:
        return surfDef

    newDef = ""
    for i, line in enumerate(lines):
        if i in longLine:
            newLine = cutLine(line, lineLength)
        else:
            newLine = line
        newDef += newLine + "\n"

    return newDef[:-1]


def cutLine(line, lineLength):
    tabNumber = 10
    while True:
        pos = line.rfind(" ")
        if pos <= lineLength:
            break

    line1 = line[0:pos]
    line2 = line[pos + 1 :]

    # second line is only spaces
    if len(line2.strip()) == 0:
        newLine = line1
    else:
        newLine = "{}\n{: <{n}}{}".format(line1, "", line2, n=tabNumber)

    return newLine
