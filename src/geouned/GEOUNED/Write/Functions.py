import math
import re

import FreeCAD

from ..Utils import Qform as Qform
from ..Utils.BasicFunctions_part1 import is_opposite, is_parallel
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
                line = f"{self.tabstr}{data_str}"
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

    # TODO check this is used
    def del_ast_char(self):
        self.str = self.str[0:-1]
        self.__leftspace__ += 1

    def wrap_line(self, offset=0):

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
                newline += f"\n{self.tabstr}"
            else:
                newline += f"\n{self.tabstr}{self.str[init:]}"
                break

        self.str = newline


def change_surf_sign(surf, Seq):
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
                change_surf_sign(surf, e)


def write_mcnp_cell_def(definition, tabspace=0, offset=0):
    sdef = CellString(tabspace=tabspace)
    str_def = remove_redundant(write_sequence_mcnp(definition))
    sdef.add(str_def)
    sdef.wrap_line(offset)
    return sdef.str


def write_serpent_cell_def(definition, tabspace=0, offset=0):
    sdef = CellString(tabspace=tabspace)
    str_def = remove_redundant(write_sequence_serpent(definition))
    sdef.add(str_def)
    sdef.wrap_line(offset)
    return sdef.str


def write_phits_cell_def(definition, tabspace=0, offset=0):
    sdef = CellString(tabspace=tabspace)
    str_def = remove_redundant(write_sequence_phits(definition))
    sdef.add(str_def)
    sdef.wrap_line(offset)
    return sdef.str


def write_openmc_region(definition, w_type="XML"):
    if w_type == "XML":
        return write_sequence_omc_xml(definition)
    if w_type == "PY":
        return write_sequence_omc_py(definition)


def write_sequence_mcnp(Seq):
    if Seq.level == 0:
        if Seq.operator == "AND":
            line = f"({' '.join(map(str, Seq.elements))})"
        else:
            line = f"({':'.join(map(str, Seq.elements))})"
    else:
        terms = []
        for e in Seq.elements:
            if type(e) is int:
                terms.append(str(e))
            else:
                terms.append(write_sequence_mcnp(e))

        if Seq.operator == "AND":
            line = f"({' '.join(terms)})"
        else:
            line = f"({':'.join(terms)})"

    return line


def write_sequence_serpent(seq):
    if seq.level == 0:
        if seq.operator == "AND":
            line = f"({' '.join(map(str, seq.elements))})"
        else:
            line = f"({':'.join(map(str, seq.elements))})"
    else:
        terms = []
        for e in seq.elements:
            if type(e) is int:
                terms.append(str(e))
            else:
                terms.append(write_sequence_mcnp(e))

        if seq.operator == "AND":
            line = f"({' '.join(terms)})"
        else:
            line = f"({':'.join(terms)})"

    return line


def write_sequence_phits(seq):
    if seq.level == 0:
        if seq.operator == "AND":
            line = f"({' '.join(map(str, seq.elements))})"
        else:
            line = f"({':'.join(map(str, seq.elements))})"
    else:
        terms = []
        for e in seq.elements:
            if type(e) is int:
                terms.append(str(e))
            else:
                terms.append(write_sequence_phits(e))

        if seq.operator == "AND":
            line = f"({' '.join(terms)})"
        else:
            line = f"({':'.join(terms)})"

    return line


def write_sequence_omc_xml(seq):
    if seq.level == 0:
        if seq.operator == "AND":
            line = f"({' '.join(map(str, seq.elements))})"
        else:
            line = f"({' | '.join(map(str, seq.elements))})"
    else:
        terms = []
        for e in seq.elements:
            if type(e) is int:
                terms.append(str(e))
            else:
                terms.append(write_sequence_omc_xml(e))

        if seq.operator == "AND":
            line = f"({' '.join(terms)})"
        else:
            line = f"({' | '.join(terms)})"
    return line


def write_sequence_omc_py(seq, prefix="S"):

    strSurf = lambda surf: (f"-{prefix}{-surf}" if surf < 0 else f"+{prefix}{surf}")

    if seq.level == 0:
        if seq.operator == "AND":
            line = f"({' & '.join(map(strSurf, seq.elements))})"
        else:
            line = f"({' | '.join(map(strSurf, seq.elements))})"
    else:
        terms = []
        for e in seq.elements:
            if type(e) is int:
                terms.append(strSurf(e))
            else:
                terms.append(write_sequence_omc_py(e))

        if seq.operator == "AND":
            line = f"({' & '.join(terms)})"
        else:
            line = f"({' | '.join(terms)})"
    return line


def mcnp_surface(id, surface_type, surf, tolerances, options, numeric_format):
    mcnp_def = ""

    if surface_type == "Plane":
        if surf.pointDef and options.prnt3PPlane:
            P1 = surf.Points[0]
            P2 = surf.Points[1]
            P3 = surf.Points[2]
            mcnp_def = """{:<6d} P   {P1[0]:{d}} {P1[1]:{d}} {P1[2]:{d}} 
           {P2[0]:{d}} {P2[1]:{d}} {P2[2]:{d}}
           {P3[0]:{d}} {P3[1]:{d}} {P3[2]:{d}}""".format(
                id, P1=P1 / 10, P2=P2 / 10, P3=P3 / 10, d=numeric_format.P_d
            )
        else:
            A = surf.Axis.x
            B = surf.Axis.y
            C = surf.Axis.z
            D = surf.Axis.dot(surf.Position)
            if surf.Axis.isEqual(FreeCAD.Vector(1, 0, 0), tolerances.pln_angle):
                mcnp_def = "{:<6d} PX  {:{x}}".format(
                    id, D / 10.0, x=numeric_format.P_xyz
                )
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 1, 0), tolerances.pln_angle):
                mcnp_def = "{:<6d} PY  {:{y}}".format(
                    id, D / 10.0, y=numeric_format.P_xyz
                )
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 0, 1), tolerances.pln_angle):
                mcnp_def = "{:<6d} PZ  {:{z}}".format(
                    id, D / 10.0, z=numeric_format.P_xyz
                )
            else:
                mcnp_def = "{:<6d} P   {:{abc}} {:{abc}} {:{abc}} {:{d}}".format(
                    id,
                    A,
                    B,
                    C,
                    D / 10.0,
                    abc=numeric_format.P_abc,
                    d=numeric_format.P_d,
                )

    elif surface_type == "Cylinder":
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        Pos = surf.Center * 0.1
        rad = surf.Radius * 0.1
        if is_parallel(Dir, FreeCAD.Vector(1, 0, 0), tolerances.angle):
            if Pos.y == 0.0 and Pos.z == 0.0:
                mcnp_def = "{:<6d} CX  {:{r}}".format(id, rad, r=numeric_format.C_r)
            else:
                mcnp_def = "{:<6d} C/X  {:{yz}} {:{yz}} {:{r}}".format(
                    id, Pos.y, Pos.z, rad, yz=numeric_format.C_xyz, r=numeric_format.C_r
                )
        elif is_parallel(Dir, FreeCAD.Vector(0, 1, 0), tolerances.angle):
            if Pos.x == 0.0 and Pos.z == 0.0:
                mcnp_def = "{:<6d} CY  {:{r}}".format(id, rad, r=numeric_format.C_r)
            else:
                mcnp_def = "{:<6d} C/Y  {:{xz}} {:{xz}} {:{r}}".format(
                    id, Pos.x, Pos.z, rad, xz=numeric_format.C_xyz, r=numeric_format.C_r
                )
        elif is_parallel(Dir, FreeCAD.Vector(0, 0, 1), tolerances.angle):
            if Pos.y == 0.0 and Pos.x == 0.0:
                mcnp_def = "{:<6d} CZ  {:{r}}".format(id, rad, r=numeric_format.C_r)
            else:
                mcnp_def = "{:<6d} C/Z  {:{xy}} {:{xy}} {:{r}}".format(
                    id, Pos.x, Pos.y, rad, xy=numeric_format.C_xyz, r=numeric_format.C_r
                )
        else:
            # Is not still working fine
            Q = Qform.q_form_cyl(Dir, Pos, rad)
            mcnp_def = """\
{:<6d} GQ  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id,
                v=Q,
                aTof=numeric_format.GQ_1to6,
                gToi=numeric_format.GQ_7to9,
                j=numeric_format.GQ_10,
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
        #    mcnp_def='%i  RCC  %13.7E %13.7E %13.7E %13.7E\n       %13.7E %13.7E %13.7E' %(id,Vx,Vy,Vz,Hx,Hy,Hz,rad)

    elif surface_type == "Cone":
        Apex = surf.Apex * 0.1
        Dir = surf.Axis * 0.1
        tan = math.tan(surf.SemiAngle)
        X_dir = FreeCAD.Vector(1, 0, 0)
        Y_dir = FreeCAD.Vector(0, 1, 0)
        Z_dir = FreeCAD.Vector(0, 0, 1)
        if is_parallel(Dir, X_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, X_dir, tolerances.angle):
                sheet = -1
            if Apex.y == 0.0 and Apex.z == 0.0:
                mcnp_def = "{:<6d} KX  {:{x}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    tan**2,
                    sheet,
                    x=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                mcnp_def = "{:<6d} K/X  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
        elif is_parallel(Dir, Y_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, Y_dir, tolerances.angle):
                sheet = -1
            if Apex.x == 0.0 and Apex.z == 0.0:
                mcnp_def = "{:<6d} KY  {:{y}} {:{t2}} {}".format(
                    id,
                    Apex.y,
                    tan**2,
                    sheet,
                    y=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                mcnp_def = "{:<6d} K/Y  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
        elif is_parallel(Dir, Z_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, Z_dir, tolerances.angle):
                sheet = -1
            if Apex.x == 0.0 and Apex.y == 0.0:
                mcnp_def = "{:<6d} KZ  {:{z}} {:{t2}} {}".format(
                    id,
                    Apex.z,
                    tan**2,
                    sheet,
                    z=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                mcnp_def = "{:<6d} K/Z  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
        else:
            Q = Qform.q_form_cone(Dir, Apex, tan)
            mcnp_def = """\
{:<6d} GQ  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id,
                v=Q,
                aTof=numeric_format.GQ_1to6,
                gToi=numeric_format.GQ_7to9,
                j=numeric_format.GQ_10,
            )

    elif surface_type == "Sphere":
        # corresponding logic
        rad = surf.Radius * 0.1
        pnt = surf.Center * 0.1
        if pnt.isEqual(FreeCAD.Vector(0, 0, 0), tolerances.sph_distance):
            mcnp_def = "{:<6d} SO  {:{r}}".format(id, rad, r=numeric_format.S_r)
        else:
            mcnp_def = "{:<6d} S  {:{xyz}} {:{xyz}} {:{xyz}} {:{r}}".format(
                id,
                pnt.x,
                pnt.y,
                pnt.z,
                rad,
                xyz=numeric_format.S_xyz,
                r=numeric_format.S_r,
            )

    elif surface_type == "Torus":
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        Pos = surf.Center * 0.1
        radMaj = surf.MajorRadius * 0.1
        radMin = surf.MinorRadius * 0.1
        if is_parallel(Dir, FreeCAD.Vector(1, 0, 0), tolerances.angle):
            mcnp_def = """\
{:<6d} TX  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id,
                Pos.x,
                Pos.y,
                Pos.z,
                radMaj,
                radMin,
                radMin,
                xyz=numeric_format.T_xyz,
                r=numeric_format.T_r,
            )
        elif is_parallel(Dir, FreeCAD.Vector(0, 1, 0), tolerances.angle):
            mcnp_def = """\
{:<6d} TY  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id,
                Pos.x,
                Pos.y,
                Pos.z,
                radMaj,
                radMin,
                radMin,
                xyz=numeric_format.T_xyz,
                r=numeric_format.T_r,
            )
        elif is_parallel(Dir, FreeCAD.Vector(0, 0, 1), tolerances.angle):
            mcnp_def = """\
{:<6d} TZ  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id,
                Pos.x,
                Pos.y,
                Pos.z,
                radMaj,
                radMin,
                radMin,
                xyz=numeric_format.T_xyz,
                r=numeric_format.T_r,
            )

    return trim(mcnp_def, 80)


def open_mc_surface(
    surface_type, surf, tolerances, numeric_format, out_xml=True, quadricForm=False
):
    if surface_type == "Plane":
        A = surf.Axis.x
        B = surf.Axis.y
        C = surf.Axis.z
        D = surf.Axis.dot(surf.Position) * 0.1
        if surf.Axis.isEqual(FreeCAD.Vector(1, 0, 0), tolerances.pln_angle):
            if out_xml:
                omc_surf = "x-plane"
                coeffs = "{:{x}}".format(D, x=numeric_format.P_xyz)
            else:
                omc_surf = "XPlane"
                coeffs = f"x0={D}"

        elif surf.Axis.isEqual(FreeCAD.Vector(0, 1, 0), tolerances.pln_angle):
            if out_xml:
                omc_surf = "y-plane"
                coeffs = "{:{x}}".format(D, x=numeric_format.P_xyz)
            else:
                omc_surf = "YPlane"
                coeffs = f"y0={D}"

        elif surf.Axis.isEqual(FreeCAD.Vector(0, 0, 1), tolerances.pln_angle):
            if out_xml:
                omc_surf = "z-plane"
                coeffs = "{:{x}}".format(D, x=numeric_format.P_xyz)
            else:
                omc_surf = "ZPlane"
                coeffs = f"z0={D}"

        else:
            if out_xml:
                omc_surf = "plane"
                coeffs = "{:{abc}} {:{abc}} {:{abc}} {:{d}}".format(
                    A, B, C, D, abc=numeric_format.P_abc, d=numeric_format.P_d
                )
            else:
                omc_surf = "Plane"
                coeffs = f"a={A},b={B},c={C},d={D}"

    elif surface_type == "Cylinder":
        pos = surf.Center * 0.1
        Rad = surf.Radius * 0.1
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()

        if is_parallel(Dir, FreeCAD.Vector(1, 0, 0), tolerances.angle):
            if out_xml:
                omc_surf = "x-cylinder"
                coeffs = "{:{xy}} {:{xy}} {:{r}}".format(
                    pos.y, pos.z, Rad, xy=numeric_format.C_xyz, r=numeric_format.C_r
                )
            else:
                omc_surf = "XCylinder"
                coeffs = f"y0={pos.y},z0={pos.z},r={Rad}"

        elif is_parallel(Dir, FreeCAD.Vector(0, 1, 0), tolerances.angle):
            if out_xml:
                omc_surf = "y-cylinder"
                coeffs = "{:{xy}} {:{xy}} {:{r}}".format(
                    pos.x, pos.z, Rad, xy=numeric_format.C_xyz, r=numeric_format.C_r
                )
            else:
                omc_surf = "YCylinder"
                coeffs = f"x0={pos.x},z0={pos.z},r={Rad}"

        elif is_parallel(Dir, FreeCAD.Vector(0, 0, 1), tolerances.angle):
            if out_xml:
                omc_surf = "z-cylinder"
                coeffs = "{:{xy}} {:{xy}} {:{r}}".format(
                    pos.x, pos.y, Rad, xy=numeric_format.C_xyz, r=numeric_format.C_r
                )
            else:
                omc_surf = "ZCylinder"
                coeffs = f"x0={pos.x},y0={pos.y},r={Rad}"

        else:
            if out_xml:
                omc_surf = "quadric"
                Q = Qform.q_form_cyl(Dir, pos, Rad)
                coeffs = "{v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}} {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}} {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}} {v[9]:{j}}".format(
                    v=Q,
                    aTof=numeric_format.GQ_1to6,
                    gToi=numeric_format.GQ_7to9,
                    j=numeric_format.GQ_10,
                )
            else:
                if quadricForm:
                    omc_surf = "Quadric"
                    Q = Qform.q_form_cyl(Dir, pos, Rad)
                    coeffs = "a={v[0]},b={v[1]},c={v[2]},d={v[3]},e={v[4]},f={v[5]},g={v[6]},h={v[7]},j={v[8]},k={v[9]}".format(
                        v=Q
                    )
                else:
                    omc_surf = "Cylinder"
                    coeffs = "x0={},y0={},z0={},r={},dx={},dy={},dz={}".format(
                        pos.x, pos.y, pos.z, Rad, Dir.x, Dir.y, Dir.z
                    )

    elif surface_type == "Cone":
        Apex = surf.Apex * 0.1
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        tan = math.tan(surf.SemiAngle)
        tan2 = tan * tan

        X_dir = FreeCAD.Vector(1, 0, 0)
        Y_dir = FreeCAD.Vector(0, 1, 0)
        Z_dir = FreeCAD.Vector(0, 0, 1)

        if is_parallel(Dir, X_dir, tolerances.angle):
            if out_xml:
                omc_surf = "x-cone"
                coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{t2}}".format(
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan2,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                omc_surf = "XCone"
                coeffs = f"x0={Apex.x},y0={Apex.y},z0={Apex.z},r2={tan2}"

        elif is_parallel(Dir, Y_dir, tolerances.angle):
            if out_xml:
                omc_surf = "y-cone"
                coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{t2}}".format(
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan2,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                omc_surf = "YCone"
                coeffs = f"x0={Apex.x},y0={Apex.y},z0={Apex.z},r2={tan2}"

        elif is_parallel(Dir, Z_dir, tolerances.angle):
            if out_xml:
                omc_surf = "z-cone"
                coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{t2}}".format(
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan2,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                omc_surf = "ZCone"
                coeffs = f"x0={Apex.x},y0={Apex.y},z0={Apex.z},r2={tan2}"

        else:
            if out_xml:
                omc_surf = "quadric"
                Q = Qform.q_form_cone(Dir, Apex, tan)
                coeffs = "{v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}} {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}} {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}} {v[9]:{j}}".format(
                    v=Q,
                    aTof=numeric_format.GQ_1to6,
                    gToi=numeric_format.GQ_7to9,
                    j=numeric_format.GQ_10,
                )
            else:
                if quadricForm:
                    omc_surf = "Quadric"
                    Q = Qform.q_form_cone(Dir, Apex, tan)
                    coeffs = "a={v[0]},b={v[1]},c={v[2]},d={v[3]},e={v[4]},f={v[5]},g={v[6]},h={v[7]},j={v[8]},k={v[9]}".format(
                        v=Q
                    )
                else:
                    omc_surf = "Cone"
                    coeffs = "x0={},y0={},z0={},r2={},dx={},dy={},dz={}".format(
                        Apex.x, Apex.y, Apex.z, tan2, Dir.x, Dir.y, Dir.z
                    )

    elif surface_type == "Sphere":
        Center = surf.Center * 0.1
        Rad = surf.Radius * 0.1
        if out_xml:
            omc_surf = "sphere"
            coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{r}}".format(
                Center.x,
                Center.y,
                Center.z,
                Rad,
                xyz=numeric_format.S_xyz,
                r=numeric_format.S_r,
            )
        else:
            omc_surf = "Sphere"
            coeffs = f"x0={Center.x},y0={Center.y},z0={Center.z},r={Rad}"

    elif surface_type == "Torus":
        Center = surf.Center * 0.1
        minRad = surf.MinorRadius * 0.1
        majRad = surf.MajorRadius * 0.1
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        if out_xml:
            coeffs = "{:{xyz}} {:{xyz}} {:{xyz}} {:{r}} {:{r}} {:{r}}".format(
                Center.x,
                Center.y,
                Center.z,
                majRad,
                minRad,
                minRad,
                xyz=numeric_format.T_xyz,
                r=numeric_format.T_r,
            )
        else:
            coeffs = "x0={},y0={},z0={},r={},r1={},r2={}".format(
                Center.x, Center.y, Center.z, majRad, minRad, minRad
            )

        if is_parallel(Dir, FreeCAD.Vector(1, 0, 0), tolerances.angle):
            omc_surf = "x-torus" if out_xml else "XTorus"
        elif is_parallel(Dir, FreeCAD.Vector(0, 1, 0), tolerances.angle):
            omc_surf = "y-torus" if out_xml else "YTorus"
        elif is_parallel(Dir, FreeCAD.Vector(0, 0, 1), tolerances.angle):
            omc_surf = "z-torus" if out_xml else "ZTorus"
        else:
            omc_surf = None

    if out_xml:
        coeffs = " ".join(coeffs.split())
    return omc_surf, coeffs


def serpent_surface(id, surface_type, surf, tolerances, numeric_format, options):
    serpent_def = ""

    if surface_type == "Plane":
        if surf.pointDef and options.prnt3PPlane:
            P1 = surf.Points[0]
            P2 = surf.Points[1]
            P3 = surf.Points[2]
            serpent_def = f"surf {id} plane {P1.x/10:{numeric_format.P_d}} {P1.y/10:{numeric_format.P_d}} {P1.z/10:{numeric_format.P_d}}\n"
            serpent_def += f"      {P2.x/10:{numeric_format.P_d}} {P2.y/10:{numeric_format.P_d}} {P2.z/10:{numeric_format.P_d}}\n"
            serpent_def += f"      {P3.x/10:{numeric_format.P_d}} {P3.y/10:{numeric_format.P_d}} {P3.z/10:{numeric_format.P_d}}"

        else:
            A = surf.Axis.x
            B = surf.Axis.y
            C = surf.Axis.z
            D = surf.Axis.dot(surf.Position)
            if surf.Axis.isEqual(FreeCAD.Vector(1, 0, 0), tolerances.pln_angle):
                serpent_def = f"surf {id} px {D/10:{numeric_format.P_xyz}}"
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 1, 0), tolerances.pln_angle):
                serpent_def = f"surf {id} py {D/10:{numeric_format.P_xyz}}"
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 0, 1), tolerances.pln_angle):
                serpent_def = f"surf {id} pz {D/10:{numeric_format.P_xyz}}"
            else:
                serpent_def = f"surf {id} plane {A:{numeric_format.P_d}} {B:{numeric_format.P_d}} {C:{numeric_format.P_d}} {D/10:{numeric_format.P_d}}"

    elif surface_type == "Cylinder":
        Dir = surf.Axis
        Dir.normalize()
        Pos = surf.Center * 0.1
        rad = surf.Radius * 0.1
        if is_parallel(Dir, FreeCAD.Vector(1, 0, 0), tolerances.angle):
            serpent_def = f"surf {id} cylx {Pos.y:{numeric_format.C_xyz}} {Pos.z:{numeric_format.C_xyz}} {rad:{numeric_format.C_r}}"
        elif is_parallel(Dir, FreeCAD.Vector(0, 1, 0), tolerances.angle):
            serpent_def = f"surf {id} cyly {Pos.x:{numeric_format.C_xyz}} {Pos.z:{numeric_format.C_xyz}} {rad:{numeric_format.C_r}}"
        elif is_parallel(Dir, FreeCAD.Vector(0, 0, 1), tolerances.angle):
            serpent_def = f"surf {id} cylz {rad:{numeric_format.C_r}}"
        else:
            # Is not still working fine
            Q = Qform.q_form_cyl(Dir, Pos, rad)
            serpent_def = """\
surf quadratic  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id,
                v=Q,
                aTof=numeric_format.GQ_1to6,
                gToi=numeric_format.GQ_7to9,
                j=numeric_format.GQ_10,
            )

    elif surface_type == "Cone":
        Apex = surf.Apex * 0.1
        Dir = surf.Axis * 0.1
        tan = math.tan(surf.SemiAngle)
        X_dir = FreeCAD.Vector(1, 0, 0)
        Y_dir = FreeCAD.Vector(0, 1, 0)
        Z_dir = FreeCAD.Vector(0, 0, 1)

        # Need to check this
        # Serpent has no specific card for cone at origin, explicit origin only

        if is_parallel(Dir, X_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, X_dir, tolerances.angle):
                sheet = -1
            serpent_def = "surf ckx {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                id,
                Apex.x,
                Apex.y,
                Apex.z,
                tan**2,
                sheet,
                xyz=numeric_format.K_xyz,
                t2=numeric_format.K_tan2,
            )
        elif is_parallel(Dir, Y_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, Y_dir, tolerances.angle):
                sheet = -1
            serpent_def = "surf cky {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                id,
                Apex.x,
                Apex.y,
                Apex.z,
                tan**2,
                sheet,
                xyz=numeric_format.K_xyz,
                t2=numeric_format.K_tan2,
            )
        elif is_parallel(Dir, Z_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, Z_dir, tolerances.angle):
                sheet = -1
            serpent_def = "surf ckz {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                id,
                Apex.x,
                Apex.y,
                Apex.z,
                tan**2,
                sheet,
                xyz=numeric_format.K_xyz,
                t2=numeric_format.K_tan2,
            )
        else:
            Q = Qform.q_form_cone(Dir, Apex, tan)

    elif surface_type == "Sphere":
        rad = surf.Radius * 0.1
        pnt = surf.Center * 0.1
        # Serpent has only explicit spheres at the origin
        serpent_def = f"surf {id} sph {pnt.x:{numeric_format.S_xyz}} {pnt.y:{numeric_format.S_xyz}} {pnt.z:{numeric_format.S_xyz}} {rad:{numeric_format.S_r}}"

    elif surface_type == "Torus":
        Dir = surf.Axis
        Dir.normalize()
        Pos = surf.Center * 0.1
        radMaj = surf.MajorRadius * 0.1
        radMin = surf.MinorRadius * 0.1
        if is_parallel(Dir, FreeCAD.Vector(1, 0, 0), tolerances.angle):
            serpent_def = f"surf {id} torx {Pos.x:{numeric_format.T_xyz}} {Pos.y:{numeric_format.T_xyz}} {Pos.z:{numeric_format.T_xyz}}\n"
            serpent_def += f"      {radMaj:{numeric_format.T_r}} {radMin:{numeric_format.T_r}} {radMin:{numeric_format.T_r}}"
        elif is_parallel(Dir, FreeCAD.Vector(0, 1, 0), tolerances.angle):
            serpent_def = f"surf {id} tory {Pos.x:{numeric_format.T_xyz}} {Pos.y:{numeric_format.T_xyz}} {Pos.z:{numeric_format.T_xyz}}\n"
            serpent_def += f"      {radMaj:{numeric_format.T_r}} {radMin:{numeric_format.T_r}} {radMin:{numeric_format.T_r}}"
        elif is_parallel(Dir, FreeCAD.Vector(0, 0, 1), tolerances.angle):
            serpent_def = f"surf {id} torz {Pos.x:{numeric_format.T_xyz}} {Pos.y:{numeric_format.T_xyz}} {Pos.z:{numeric_format.T_xyz}}\n"
            serpent_def += f"      {radMaj:{numeric_format.T_r}} {radMin:{numeric_format.T_r}} {radMin:{numeric_format.T_r}}"

    return serpent_def


def phits_surface(id, surface_type, surf, tolerances, numeric_format, options):
    phits_def = ""

    if surface_type == "Plane":
        if surf.pointDef and options.prnt3PPlane:
            P1 = surf.Points[0]
            P2 = surf.Points[1]
            P3 = surf.Points[2]
            phits_def = """{:<6d} P   {P1[0]:{d}} {P1[1]:{d}} {P1[2]:{d}} 
           {P2[0]:{d}} {P2[1]:{d}} {P2[2]:{d}}
           {P3[0]:{d}} {P3[1]:{d}} {P3[2]:{d}}""".format(
                id, P1=P1 / 10, P2=P2 / 10, P3=P3 / 10, d=numeric_format.P_d
            )
        else:
            A = surf.Axis.x
            B = surf.Axis.y
            C = surf.Axis.z
            D = surf.Axis.dot(surf.Position)
            if surf.Axis.isEqual(FreeCAD.Vector(1, 0, 0), tolerances.pln_angle):
                phits_def = "{:<6d} PX  {:{x}}".format(
                    id, D / 10.0, x=numeric_format.P_xyz
                )
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 1, 0), tolerances.pln_angle):
                phits_def = "{:<6d} PY  {:{y}}".format(
                    id, D / 10.0, y=numeric_format.P_xyz
                )
            elif surf.Axis.isEqual(FreeCAD.Vector(0, 0, 1), tolerances.pln_angle):
                phits_def = "{:<6d} PZ  {:{z}}".format(
                    id, D / 10.0, z=numeric_format.P_xyz
                )
            else:
                phits_def = "{:<6d} P   {:{abc}} {:{abc}} {:{abc}} {:{d}}".format(
                    id,
                    A,
                    B,
                    C,
                    D / 10.0,
                    abc=numeric_format.P_abc,
                    d=numeric_format.P_d,
                )

    elif surface_type == "Cylinder":
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        Pos = surf.Center * 0.1
        rad = surf.Radius * 0.1
        if is_parallel(Dir, FreeCAD.Vector(1, 0, 0), tolerances.angle):
            if Pos.y == 0.0 and Pos.z == 0.0:
                phits_def = "{:<6d} CX  {:{r}}".format(id, rad, r=numeric_format.C_r)
            else:
                phits_def = "{:<6d} C/X  {:{yz}} {:{yz}} {:{r}}".format(
                    id, Pos.y, Pos.z, rad, yz=numeric_format.C_xyz, r=numeric_format.C_r
                )
        elif is_parallel(Dir, FreeCAD.Vector(0, 1, 0), tolerances.angle):
            if Pos.x == 0.0 and Pos.z == 0.0:
                phits_def = "{:<6d} CY  {:{r}}".format(id, rad, r=numeric_format.C_r)
            else:
                phits_def = "{:<6d} C/Y  {:{xz}} {:{xz}} {:{r}}".format(
                    id, Pos.x, Pos.z, rad, xz=numeric_format.C_xyz, r=numeric_format.C_r
                )
        elif is_parallel(Dir, FreeCAD.Vector(0, 0, 1), tolerances.angle):
            if Pos.y == 0.0 and Pos.x == 0.0:
                phits_def = "{:<6d} CZ  {:{r}}".format(id, rad, r=numeric_format.C_r)
            else:
                phits_def = "{:<6d} C/Z  {:{xy}} {:{xy}} {:{r}}".format(
                    id, Pos.x, Pos.y, rad, xy=numeric_format.C_xyz, r=numeric_format.C_r
                )
        else:
            # Is not still working fine
            Q = Qform.q_form_cyl(Dir, Pos, rad)
            phits_def = """\
{:<6d} GQ  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id,
                v=Q,
                aTof=numeric_format.GQ_1to6,
                gToi=numeric_format.GQ_7to9,
                j=numeric_format.GQ_10,
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
        #    mcnp_def='%i  RCC  %13.7E %13.7E %13.7E %13.7E\n       %13.7E %13.7E %13.7E' %(id,Vx,Vy,Vz,Hx,Hy,Hz,rad)

    elif surface_type == "Cone":
        Apex = surf.Apex * 0.1
        Dir = surf.Axis * 0.1
        tan = math.tan(surf.SemiAngle)
        X_dir = FreeCAD.Vector(1, 0, 0)
        Y_dir = FreeCAD.Vector(0, 1, 0)
        Z_dir = FreeCAD.Vector(0, 0, 1)
        if is_parallel(Dir, X_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, X_dir, tolerances.angle):
                sheet = -1
            if Apex.y == 0.0 and Apex.z == 0.0:
                phits_def = "{:<6d} KX  {:{x}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    tan**2,
                    sheet,
                    x=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                phits_def = "{:<6d} K/X  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
        elif is_parallel(Dir, Y_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, Y_dir, tolerances.angle):
                sheet = -1
            if Apex.x == 0.0 and Apex.z == 0.0:
                phits_def = "{:<6d} KY  {:{y}} {:{t2}} {}".format(
                    id,
                    Apex.y,
                    tan**2,
                    sheet,
                    y=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                phits_def = "{:<6d} K/Y  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
        elif is_parallel(Dir, Z_dir, tolerances.angle):
            sheet = 1
            if is_opposite(Dir, Z_dir, tolerances.angle):
                sheet = -1
            if Apex.x == 0.0 and Apex.y == 0.0:
                phits_def = "{:<6d} KZ  {:{z}} {:{t2}} {}".format(
                    id,
                    Apex.z,
                    tan**2,
                    sheet,
                    z=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
            else:
                phits_def = "{:<6d} K/Z  {:{xyz}} {:{xyz}} {:{xyz}} {:{t2}} {}".format(
                    id,
                    Apex.x,
                    Apex.y,
                    Apex.z,
                    tan**2,
                    sheet,
                    xyz=numeric_format.K_xyz,
                    t2=numeric_format.K_tan2,
                )
        else:
            Q = Qform.q_form_cone(Dir, Apex, tan)
            phits_def = """\
{:<6d} GQ  {v[0]:{aTof}} {v[1]:{aTof}} {v[2]:{aTof}}
          {v[3]:{aTof}} {v[4]:{aTof}} {v[5]:{aTof}}
          {v[6]:{gToi}} {v[7]:{gToi}} {v[8]:{gToi}}
          {v[9]:{j}} """.format(
                id,
                v=Q,
                aTof=numeric_format.GQ_1to6,
                gToi=numeric_format.GQ_7to9,
                j=numeric_format.GQ_10,
            )

    elif surface_type == "Sphere":
        # corresponding logic
        rad = surf.Radius * 0.1
        pnt = surf.Center * 0.1
        if pnt.isEqual(FreeCAD.Vector(0, 0, 0), tolerances.sph_distance):
            phits_def = "{:<6d} SO  {:{r}}".format(id, rad, r=numeric_format.S_r)
        else:
            phits_def = "{:<6d} S  {:{xyz}} {:{xyz}} {:{xyz}} {:{r}}".format(
                id,
                pnt.x,
                pnt.y,
                pnt.z,
                rad,
                xyz=numeric_format.S_xyz,
                r=numeric_format.S_r,
            )

    elif surface_type == "Torus":
        Dir = FreeCAD.Vector(surf.Axis)
        Dir.normalize()
        Pos = surf.Center * 0.1
        radMaj = surf.MajorRadius * 0.1
        radMin = surf.MinorRadius * 0.1
        if is_parallel(Dir, FreeCAD.Vector(1, 0, 0), tolerances.angle):
            phits_def = """\
{:<6d} TX  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id,
                Pos.x,
                Pos.y,
                Pos.z,
                radMaj,
                radMin,
                radMin,
                xyz=numeric_format.T_xyz,
                r=numeric_format.T_r,
            )
        elif is_parallel(Dir, FreeCAD.Vector(0, 1, 0), tolerances.angle):
            phits_def = """\
{:<6d} TY  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id,
                Pos.x,
                Pos.y,
                Pos.z,
                radMaj,
                radMin,
                radMin,
                xyz=numeric_format.T_xyz,
                r=numeric_format.T_r,
            )
        elif is_parallel(Dir, FreeCAD.Vector(0, 0, 1), tolerances.angle):
            phits_def = """\
{:<6d} TZ  {:{xyz}} {:{xyz}} {:{xyz}}
           {:{r}} {:{r}} {:{r}}""".format(
                id,
                Pos.x,
                Pos.y,
                Pos.z,
                radMaj,
                radMin,
                radMin,
                xyz=numeric_format.T_xyz,
                r=numeric_format.T_r,
            )

    return trim(phits_def, 80)


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
            newLine = cut_line(line, lineLength)
        else:
            newLine = line
        newDef += newLine + "\n"

    return newDef[:-1]


def cut_line(line, lineLength):
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
