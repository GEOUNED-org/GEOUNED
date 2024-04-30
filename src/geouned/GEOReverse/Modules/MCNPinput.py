import math
import os
import re

import FreeCAD
import numpy as np
from numpy import linalg as LA

from .Objects import (
    CadCell,
    Plane,
    Sphere,
    Cylinder,
    Cone,
    EllipticCone,
    Hyperboloid,
    Ellipsoid,
    EllipticCylinder,
    HyperbolicCylinder,
    Paraboloid,
    Torus,
    Box,
)
from .Parser import parser as mp
from .remh import CellCardString, remove_hash


class McnpInput:
    def __init__(self, name):
        if not os.path.isfile(name):
            raise FileNotFoundError(f"File {name} does not exist")
        self.__inputcards__ = list(mp.get_cards(name))
        self.Transformations = self.__getTransList__()
        return

    def GetFilteredCells(self, Surfaces, config):
        levels, contLevels, Universes = self.GetLevelStructure()

        FilteredCells = {}

        Ustart = config["Ustart"]
        subUniverses = getSubUniverses(Ustart, Universes)
        subUniverses.add(Ustart)

        if config["levelMax"] == "all":
            levelMax = len(levels)
        else:
            levelMax = config["levelMax"] + 1

        levelUniverse = set()
        for lev in range(0, levelMax):
            for U in levels[lev]:
                levelUniverse.add(U)
        subUniverses = subUniverses.intersection(levelUniverse)

        for U in list(Universes.keys()):
            if U not in subUniverses:
                del Universes[U]

        for U in Universes.keys():
            FilteredCells[U] = selectCells(Universes[U], config)
            processSurfaces(FilteredCells[U], Surfaces)

        # change the surface name in surface dict
        newSurfaces = {}
        for k in Surfaces.keys():
            newkey = Surfaces[k].id
            newSurfaces[newkey] = Surfaces[k]

        for U, universe in FilteredCells.items():
            substituteLikeCell(universe, newSurfaces)

            # set cell as CAD cell Object
            for cname, c in universe.items():
                # print(cname,c.geom.str)
                universe[cname] = CadCell(c)

        return levels, FilteredCells, newSurfaces

    def GetLevelStructure(self):
        containers = []
        Universe_dict = {}
        cellCards = {}

        for c in self.__inputcards__:
            if c.ctype != mp.CID.cell:
                continue
            c.get_values()
            cstr = CellCardString("".join(c.lines))

            if cstr.TRCL:
                cstr.TRCL = TransformationMatrix(cstr.TRCL, self.Transformations)
            if cstr.TR:
                cstr.TR = TransformationMatrix(cstr.TR, self.Transformations)
            cellCards[c.name] = cstr

        setExplicitCellDefinition(cellCards)
        for cname, c in cellCards.items():
            if c.U is None:
                c.U = 0
            if c.U not in Universe_dict.keys():
                Universe_dict[c.U] = {}
            Universe_dict[c.U].update({c.name: c})

            if c.FILL:
                containers.append(c)

        currentLevel = [0]
        nextLevel = []
        contLevel = {0: [(0, 0)]}
        univLevel = {0: {0}}
        level = 0

        while True:
            level += 1
            contLevel[level] = []
            univLevel[level] = set()
            for c in reversed(containers):
                if c.U in currentLevel:
                    c.Level = level
                    nextLevel.append(c.FILL)
                    contLevel[level].append((c.U, c.name))
                    univLevel[level].add(c.FILL)

                    containers.remove(c)

            if nextLevel == []:
                break
            currentLevel = nextLevel
            nextLevel = []

        return univLevel, contLevel, Universe_dict

    def GetCells(self, U=None, Fill=None):
        cell_cards = {}
        for c in self.__inputcards__:
            if c.ctype != mp.CID.cell:
                continue
            c.get_values()
            U_cell = c.get_u()
            Fill_cell = c.get_f()
            if U is None and Fill is None:
                cell_cards[c.name] = c
            elif U_cell == U and U is not None:
                cell_cards[c.name] = c
            elif Fill_cell == Fill and Fill is not None:
                cell_cards[c.name] = c

        return cell_cards

    def GetSurfaces(self, scale=1.0):
        surf_cards = {}
        number = 1
        for c in self.__inputcards__:
            if c.ctype != mp.CID.surface:
                continue
            c.get_values()
            c.TR = TransformationMatrix(c.TR, self.Transformations)
            surf_cards[c.name] = (c.stype, c.scoefs, c.TR, number)
            number += 1

        # return surface as surface Objects type
        return Get_primitive_surfaces(surf_cards, scale)

    def __getTransList__(self):
        trl = {}
        for c in self.__inputcards__:
            if c.ctype != mp.CID.data:
                continue
            c.get_values()
            if c.dtype == "TRn":
                trValues = []
                for v in c.values[1:]:
                    trValues.append(v[0])
                trl[c.name] = getTransMatrix(trValues, c.unit)
        return trl


# fmt: off
def getTransMatrix(trsf, unit="", scale=10.0):

    if len(trsf) == 3:
        trsfMat = FreeCAD.Matrix(
            1, 0, 0, trsf[0] * scale,
            0, 1, 0, trsf[1] * scale,
            0, 0, 1, trsf[2] * scale,
            0, 0, 0, 1,
        )
    else:
        if unit == "*":
            coeff = tuple(map(math.radians, trsf[3:12]))
            coeff = tuple(map(math.cos, coeff))
        else:
            coeff = trsf[3:12]

        trsfMat = FreeCAD.Matrix(
            coeff[0], coeff[3], coeff[6], trsf[0] * scale,
            coeff[1], coeff[4], coeff[7], trsf[1] * scale,
            coeff[2], coeff[5], coeff[8], trsf[2] * scale,
            0, 0, 0, 1,
        )
    return trsfMat    
# fmt: on


def substituteLikeCell(universe, Surfaces):
    number = re.compile(r"\#?\s*\d+")
    newId = len(Surfaces)

    # create new cells objects for like cell cells
    for c in universe.values():
        if c.likeCell:
            c.geom = universe[c.likeCell].geom.copy()
        if not c.TRCL:
            continue  # likebut cell should have TRCL card

    # transform change cell the parameters if needed
    for c in universe.values():
        if not c.TRCL:
            continue
        cellSurf = c.geom.getSurfacesNumbers()
        surfDict = {}
        for surf in cellSurf:
            newId += 1
            surfDict[surf] = newId
            Surfaces[newId] = Surfaces[surf].copy()
            Surfaces[newId].id = newId
            Surfaces[newId].transform(c.TRCL)
        if c.FILL:
            c.TR = c.TRCL
        c.TRCL = None

        # rename the surface of the like cell
        pos = 0
        while True:
            m = number.search(c.geom.str, pos)
            if not m:
                break
            if "#" in m.group():
                pos = m.end()
                continue
            surf = int(m.group())
            pos = c.geom.replace(surf, surfDict[surf], pos)
            if pos < 0:
                break


def selectCells(cellList, config):
    selected = {}
    # options are 'all' material
    if config["mat"][0] == "all":
        if config["cell"][0] == "all":
            selected = cellList
        elif config["cell"][0] == "exclude":
            for name, c in cellList.items():
                if name not in config["cell"][1]:
                    selected[name] = c
        elif config["cell"][0] == "include":
            for name, c in cellList.items():
                if name in config["cell"][1]:
                    selected[name] = c

    # options are 'exclude' material
    elif config["mat"][0] == "exclude":
        if config["cell"][0] == "all":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT not in config["mat"][1]:
                        selected[name] = c
                else:
                    selected[name] = (
                        c  # Fill cell are not tested against material number
                    )
        elif config["cell"][0] == "exclude":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT not in config["mat"][1]:
                        if name not in config["cell"][1]:
                            selected[name] = c
                else:
                    if name not in config["cell"][1]:
                        selected[name] = (
                            c  # Fill cell are not tested against material number
                        )
        elif config["cell"][0] == "include":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT not in config["mat"][1]:
                        if name in config["cell"][1]:
                            selected[name] = c
                else:
                    if name in config["cell"][1]:
                        selected[name] = (
                            c  # Fill cell are not tested against material number
                        )

    # options are 'include' material
    elif config["mat"][0] == "include":
        if config["cell"][0] == "all":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT in config["mat"][1]:
                        selected[name] = c
                else:
                    selected[name] = (
                        c  # Fill cell are not tested against material number
                    )
        elif config["cell"][0] == "exclude":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT in config["mat"][1]:
                        if name not in config["cell"][1]:
                            selected[name] = c
                else:
                    if name not in config["cell"][1]:
                        selected[name] = (
                            c  # Fill cell are not tested against material number
                        )
        elif config["cell"][0] == "include":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT in config["mat"][1]:
                        if name in config["cell"][1]:
                            selected[name] = c
                else:
                    if name in config["cell"][1]:
                        selected[name] = (
                            c  # Fill cell are not tested against material number
                        )

    # remove complementary in cell of the universe
    for cname, c in selected.items():
        c.geom = remove_hash(cellList, cname)

    if not selected:
        raise ValueError(
            "No cells selected. Check input or selection criteria in config file."
        )

    return selected


# Change implicit cell definition (like-but type or cell with #)
# to explicit cell definition (only surface number)
def setExplicitCellDefinition(universeCells):
    for cname, c in universeCells.items():
        if c.likeCell:
            lkc = universeCells[c.likeCell]
            c.geom = lkc.geom
            if not c.U:
                c.U = lkc.U
            if not c.FILL:
                c.FILL = lkc.FILL
            if not c.MAT:
                c.MAT = lkc.MAT
            if not c.TR:
                c.TR = lkc.TR
            if not c.TRCL:
                c.TRCL = lkc.TRCL
    return


def processSurfaces(UCells, Surfaces):
    number = re.compile(r"\#?\s*\d+")

    for cname, c in UCells.items():
        c.geom.remove_comments(full=True)
        pos = 0
        while True:
            m = number.search(c.geom.str, pos)
            if not m:
                break
            if "#" in m.group():
                pos = m.end()
                continue
            surf = int(m.group())
            if surf == 0:
                print(c.name)
                print(m)
                print(c.geom.str)
            pos = c.geom.replace(surf, Surfaces[surf].id, pos)


def getTransMatrix(trsf, unit="", scale=10.0):

    if len(trsf) == 3:
        trsfMat = FreeCAD.Matrix(
            1,
            0,
            0,
            trsf[0] * scale,
            0,
            1,
            0,
            trsf[1] * scale,
            0,
            0,
            1,
            trsf[2] * scale,
            0,
            0,
            0,
            1,
        )
    else:
        if unit == "*":
            coeff = tuple(map(math.radians, trsf[3:12]))
            coeff = tuple(map(math.cos, coeff))
        else:
            coeff = trsf[3:12]

        trsfMat = FreeCAD.Matrix(
            coeff[0],
            coeff[3],
            coeff[6],
            trsf[0] * scale,
            coeff[1],
            coeff[4],
            coeff[7],
            trsf[1] * scale,
            coeff[2],
            coeff[5],
            coeff[8],
            trsf[2] * scale,
            0,
            0,
            0,
            1,
        )
    return trsfMat


def TransformationMatrix(TRSF, Transformations):

    if TRSF:
        if type(TRSF) is int:
            return Transformations[TRSF]

        else:
            return getTransMatrix(TRSF)
    else:
        return TRSF


def getSubUniverses(Ustart, Universes):
    Uid = set()
    for c in Universes[Ustart].values():
        if c.FILL:
            Uid.add(c.FILL)

    AllU = Uid.copy()
    for U in Uid:
        AllU = AllU.union(getSubUniverses(U, Universes))

    return AllU


# traduce mcnp surface definition for Solid_Cell class
#  planes:
#     Stype = 'plane'
#     params = [ax,ay,az,d]
#
#  spheres:
#     Stype = 'shpere'
#     params = [cx,cy,cz,R]
#
#  cylinders:
#     Stype = 'cylinder'
#     params = [[px,py,pz],[vx,vy,vz],R]
#
#  cones:
#     Stype = 'cone'
#     params = [[px,py,pz],[vx,vy,vz],t,sht]
#
#  torus:
#     Stype = 'torus'
#     params = [[px,py,pz],[vx,vy,vz],ra,r]


# Return a diccionary with the corresponding surface Object
def Get_primitive_surfaces(mcnp_surfaces, scale=10.0):

    X_vec = FreeCAD.Vector(1.0, 0.0, 0.0)
    Y_vec = FreeCAD.Vector(0.0, 1.0, 0.0)
    Z_vec = FreeCAD.Vector(0.0, 0.0, 1.0)
    negX_vec = -X_vec
    negY_vec = -Y_vec
    negZ_vec = -Z_vec
    origin = FreeCAD.Vector(0.0, 0.0, 0.0)

    surfaces = {}
    for Sid in mcnp_surfaces.keys():
        MCNPtype = mcnp_surfaces[Sid][0].upper()
        MCNPparams = mcnp_surfaces[Sid][1]
        trsf = mcnp_surfaces[Sid][2]
        number = mcnp_surfaces[Sid][3]

        params = []
        Stype = None
        if MCNPtype in ("P", "PX", "PY", "PZ"):
            Stype = "plane"
            if MCNPtype == "P":
                if len(MCNPparams) == 4:
                    normal = FreeCAD.Vector(MCNPparams[0:3])
                    params = (normal, MCNPparams[3] * scale)
                else:
                    coeffs = pointsToCoeffs(MCNPparams[0:9])
                    normal = FreeCAD.Vector(coeffs[0:3])
                    point = coeffs[3] / normal.Length
                    normal.normalize()
                    params = (normal, point * scale)
            elif MCNPtype == "PX":
                params = (X_vec, MCNPparams[0] * scale)
            elif MCNPtype == "PY":
                params = (Y_vec, MCNPparams[0] * scale)
            elif MCNPtype == "PZ":
                params = (Z_vec, MCNPparams[0] * scale)

        elif MCNPtype in ["S", "SX", "SY", "SZ", "SO", "SPH"]:
            Stype = "sphere"
            if MCNPtype in ["S", "SPH"]:
                params = (
                    FreeCAD.Vector(MCNPparams[0:3]) * scale,
                    MCNPparams[3] * scale,
                )
            elif MCNPtype == "SX":
                params = (
                    FreeCAD.Vector(MCNPparams[0] * scale, 0.0, 0.0),
                    MCNPparams[1] * scale,
                )
            elif MCNPtype == "SY":
                params = (
                    FreeCAD.Vector(0.0, MCNPparams[0] * scale, 0.0),
                    MCNPparams[1] * scale,
                )
            elif MCNPtype == "SZ":
                params = (
                    FreeCAD.Vector(0.0, 0.0, MCNPparams[0] * scale),
                    MCNPparams[1] * scale,
                )
            elif MCNPtype == "SO":
                params = (origin, MCNPparams[0] * scale)

        elif MCNPtype in ["CX", "C/X", "CY", "C/Y", "CZ", "C/Z"]:
            if MCNPtype in ["CX", "CY", "CZ"]:
                R = MCNPparams[0]
            else:
                R = MCNPparams[2]
                x1 = MCNPparams[0]
                x2 = MCNPparams[1]

            Stype = "cylinder"

            if MCNPtype == "CX":
                v = X_vec
                p = origin
            elif MCNPtype == "CY":
                v = Y_vec
                p = origin
            elif MCNPtype == "CZ":
                v = Z_vec
                p = origin
            elif MCNPtype == "C/X":
                v = X_vec
                p = FreeCAD.Vector(0.0, x1, x2)
            elif MCNPtype == "C/Y":
                v = Y_vec
                p = FreeCAD.Vector(x1, 0.0, x2)
            elif MCNPtype == "C/Z":
                v = Z_vec
                p = FreeCAD.Vector(x1, x2, 0.0)

            if scale != 1.0:
                p = p.multiply(scale)
                R *= scale

            params = (p, v, R)

        elif MCNPtype in ["KX", "KY", "KZ"]:
            Stype = "cone"
            x1, t2 = MCNPparams[0:2]
            t = math.sqrt(t2)
            dblsht = True

            if len(MCNPparams) == 3:
                dblsht = False
                sht = MCNPparams[2]
                if sht == 0.0:
                    dblsht = True

            if MCNPtype == "KX":
                p = FreeCAD.Vector(x1, 0.0, 0.0)
                v = X_vec
                if not dblsht:
                    if sht < 0:
                        v = negX_vec
            elif MCNPtype == "KY":
                p = FreeCAD.Vector(0.0, x1, 0.0)
                v = Y_vec
                if not dblsht:
                    if sht < 0:
                        v = negY_vec
            elif MCNPtype == "KZ":
                p = FreeCAD.Vector(0.0, 0.0, x1)
                v = Z_vec
                if not dblsht:
                    if sht < 0:
                        v = negZ_vec

            p = p.multiply(scale)
            params = (p, v, t, dblsht)

        elif MCNPtype in ["K/X", "K/Y", "K/Z"]:
            Stype = "cone"
            p = FreeCAD.Vector(MCNPparams[0:3])
            t2 = MCNPparams[3]
            t = math.sqrt(t2)
            dblsht = True

            if len(MCNPparams) == 5:
                dblsht = False
                sht = MCNPparams[4]
                if sht == 0.0:
                    dblsht = True

            if MCNPtype == "K/X":
                v = X_vec
                if not dblsht:
                    if sht < 0:
                        v = negX_vec
            elif MCNPtype == "K/Y":
                v = Y_vec
                if not dblsht:
                    if sht < 0:
                        v = negY_vec
            elif MCNPtype == "K/Z":
                v = Z_vec
                if not dblsht:
                    if sht < 0:
                        v = negZ_vec

            p = p.multiply(scale)
            params = (p, v, t, dblsht)

        elif MCNPtype in ["TX", "TY", "TZ"]:
            Stype = "torus"
            p = FreeCAD.Vector(MCNPparams[0:3])
            Ra, Rb, Rc = MCNPparams[3:6]

            if MCNPtype == "TX":
                v = X_vec
            elif MCNPtype == "TY":
                v = Y_vec
            elif MCNPtype == "TZ":
                v = Z_vec

            if scale != 1.0:
                Ra *= scale
                Rb *= scale
                Rc *= scale
                p = p.multiply(scale)

            params = (p, v, Ra, Rb, Rc)

        elif MCNPtype == "GQ" or MCNPtype == "SQ":
            Qparams = tuple(MCNPparams[0:10])

            if MCNPtype == "SQ":
                Stype, quadric = sq2params(Qparams)
            else:
                Stype, quadric = gq2params(Qparams)

            if Stype == "cylinder":
                # p = FreeCAD.Vector(quadric[0:3])
                # v = FreeCAD.Vector(quadric[3:6])
                # R = quadric[6]
                p, v, R = quadric
                if scale != 1.0:
                    R *= scale
                    p = p.multiply(scale)

                params = (p, v, R)

            elif Stype == "cylinder_elliptic":
                p, v, radii, raxes = quadric
                if scale != 1.0:
                    radii[0] *= scale
                    radii[1] *= scale
                    p = p.multiply(scale)
                params = (p, v, radii, raxes)

            elif Stype == "cylinder_hyperbolic":
                p, v, radii, raxes = quadric
                if scale != 1.0:
                    radii[0] *= scale
                    radii[1] *= scale
                    p = p.multiply(scale)
                params = (p, v, radii, raxes)

            elif Stype == "cone":
                # p = FreeCAD.Vector(quadric[0:3])
                # v = FreeCAD.Vector(quadric[3:6])
                # t = quadric[6]
                # dblsht = quadric[7]
                p, v, t, dblsht = quadric
                if scale != 1.0:
                    p = p.multiply(scale)
                params = (p, v, t, dblsht)

            elif Stype == "cone_elliptic":
                p, v, Ra, radii, raxes, dblsht = quadric
                if scale != 1.0:
                    Ra *= scale
                    radii[0] *= scale
                    radii[1] *= scale
                    p = p.multiply(scale)
                params = (p, v, Ra, radii, raxes, dblsht)

            elif Stype == "hyperboloid":
                p, v, radii, raxes, onesht = quadric
                if scale != 1.0:
                    radii[0] *= scale
                    radii[1] *= scale
                    p = p.multiply(scale)
                params = (p, v, radii, raxes, onesht)

            elif Stype == "ellipsoid":
                p, v, radii, raxes = quadric
                if scale != 1.0:
                    radii[0] *= scale
                    radii[1] *= scale
                    p = p.multiply(scale)
                params = (p, v, radii, raxes)

            elif Stype == "paraboloid":
                p, v, focal = quadric
                if scale != 1.0:
                    focal *= scale
                    p = p.multiply(scale)

                params = (p, v, focal)

            else:
                print(Stype)
                params = None
        #                get_quadric_surface(params)

        elif MCNPtype == "X":
            if len(MCNPparams) == 2:
                Stype = "plane"
                params = (X_vec, MCNPparams[0] * scale)
            elif len(MCNPparams) == 4:
                if (abs(MCNPparams[1] - MCNPparams[3])) > 1.0e-12:
                    Stype = "cone"
                    dblsht = False
                    t = (MCNPparams[3] - MCNPparams[1]) / (
                        MCNPparams[2] - MCNPparams[0]
                    )
                    x = MCNPparams[0] - MCNPparams[1] / t
                    if (MCNPparams[0] - x) * (MCNPparams[2] - x) > 0:
                        p = FreeCAD.Vector(x, 0.0, 0.0)
                        if (MCNPparams[0] - x) > 0:
                            v = X_vec
                        else:
                            v = negX_vec
                        if scale != 1.0:
                            p *= scale
                        params = (p, v, abs(t), dblsht)
                elif abs(MCNPparams[1]) < 1.0e-12:
                    Stype = "plane"
                    if (abs(MCNPparams[0] - MCNPparams[2])) < 1.0e-12:
                        params = (X_vec, MCNPparams[0] * scale)
                else:
                    Stype = "cylinder"
                    if scale != 1.0:
                        p = p.multiply(scale)
                        R *= scale
                    params = (origin, X_vec, MCNPparams[1])
            else:
                print(
                    "not implemented surfaces defined by point with more than 2couples of value"
                )

        elif MCNPtype == "Y":
            if len(MCNPparams) == 2:
                Stype = "plane"
                params = (Y_vec, MCNPparams[0] * scale)
            elif len(MCNPparams) == 4:
                if (abs(MCNPparams[1] - MCNPparams[3])) > 1.0e-12:
                    Stype = "cone"
                    dblsht = False
                    t = (MCNPparams[3] - MCNPparams[1]) / (
                        MCNPparams[2] - MCNPparams[0]
                    )
                    y = MCNPparams[0] - MCNPparams[1] / t
                    if (MCNPparams[0] - y) * (MCNPparams[2] - y) > 0:
                        p = FreeCAD.Vector(0.0, y, 0.0)
                        if (MCNPparams[0] - y) > 0:
                            v = Y_vec
                        else:
                            v = negY_vec
                        if scale != 1.0:
                            p = p.multiply(scale)
                        params = (p, v, abs(t), dblsht)
                elif abs(MCNPparams[1]) < 1.0e-12:
                    Stype = "plane"
                    if (abs(MCNPparams[0] - MCNPparams[2])) < 1.0e-12:
                        params = (Y_vec, MCNPparams[0] * scale)
                else:
                    Stype = "cylinder"
                    if scale != 1.0:
                        p = p.multiply(scale)
                        R *= scale
                    params = (origin, Y_vec, MCNPparams[1])
            else:
                print(
                    "not implemented surfaces defined by point with more than 2couples of value"
                )

        elif MCNPtype == "Z":
            if len(MCNPparams) == 2:
                Stype = "plane"
                params = (Z_vec, MCNPparams[0] * scale)
            elif len(MCNPparams) == 4:
                if (abs(MCNPparams[1] - MCNPparams[3])) > 1.0e-12:
                    Stype = "cone"
                    dblsht = False
                    t = (MCNPparams[3] - MCNPparams[1]) / (
                        MCNPparams[2] - MCNPparams[0]
                    )
                    z = MCNPparams[0] - MCNPparams[1] / t
                    if (MCNPparams[0] - z) * (MCNPparams[2] - z) > 0:
                        p = FreeCAD.Vector(0.0, 0.0, z)
                        if (MCNPparams[0] - z) > 0:
                            v = Z_vec
                        else:
                            v = negZ_vec
                        if scale != 1.0:
                            p = p.multiply(scale)
                        params = (p, v, abs(t), dblsht)
                elif abs(MCNPparams[1]) < 1.0e-12:
                    Stype = "plane"
                    if (abs(MCNPparams[0] - MCNPparams[2])) < 1.0e-12:
                        params = (Z_vec, MCNPparams[0] * scale)
                else:
                    Stype = "cylinder"
                    if scale != 1.0:
                        p = p.multiply(scale)
                        R *= scale
                    params = (origin, Z_vec, MCNPparams[1])
            else:
                print(
                    "not implemented surfaces defined by point with more than 2couples of value"
                )

        elif MCNPtype == "BOX":
            Stype = "box"
            p = FreeCAD.Vector(MCNPparams[0:3])
            v1 = FreeCAD.Vector(MCNPparams[3:6])
            v2 = FreeCAD.Vector(MCNPparams[6:9])
            v3 = FreeCAD.Vector(MCNPparams[9:12])
            if scale != 1.0:
                p = p.multiply(scale)
                v1 = v1.multiply(scale)
                v2 = v2.multiply(scale)
                v3 = v3.multiply(scale)
            params = (p, v1, v2, v3)

        elif MCNPtype == "RPP":
            Stype = "box"
            xmin, xmax, ymin, ymax, zmin, zmax = MCNPparams[0:6]
            lx = xmax - xmin
            ly = ymax - ymin
            lz = zmax - zmin
            p = FreeCAD.Vector(xmin, ymin, zmin)
            v1 = FreeCAD.Vector(lx, 0, 0)
            v2 = FreeCAD.Vector(0, ly, 0)
            v3 = FreeCAD.Vector(0, 0, lz)
            if scale != 1.0:
                p = p.multiply(scale)
                v1 = v1.multiply(scale)
                v2 = v2.multiply(scale)
                v3 = v3.multiply(scale)
            params = (p, v1, v2, v3)

        elif MCNPtype == "RCC":
            Stype = "can"
            p = FreeCAD.Vector(MCNPparams[0:3])
            v = FreeCAD.Vector(MCNPparams[3:6])
            R = MCNPparams[6]
            if scale != 1.0:
                p = p.multiply(scale)
                v = v.multiply(scale)
                R *= scale
            params = (p, v, R)

        elif MCNPtype == "REC":
            Stype = "ecan"
            p = FreeCAD.Vector(MCNPparams[0:3])
            v = FreeCAD.Vector(MCNPparams[3:6])
            majAxis = FreeCAD.Vector(MCNPparams[6:9])
            majRad = majAxis.Length
            majAxis.normalize()

            if len(MCNPparams) == 12:
                minAxis = FreeCAD.Vector(MCNPparams[9:12])
                minRad = minAxis.Length
                minAxis.normalize()
            else:
                minRad = MCNPparams[9]
                minAxis = v.cross(majAxis) / v.Length

            if scale != 1.0:
                p = p.multiply(scale)
                v = v.multiply(scale)
                majRad *= scale
                minRad *= scale

            if minRad > majRad:
                minRad, majRad = majRad, minRad
                minAxis, majAxis = majAxis, minAxis

            params = (p, v, [minRad, majRad], [minAxis, majAxis])

        elif MCNPtype == "TRC":
            Stype = "tcone"
            p = FreeCAD.Vector(MCNPparams[0:3])
            v = FreeCAD.Vector(MCNPparams[3:6])
            R1, R2 = MCNPparams[6:8]
            if scale != 1.0:
                p = p.multiply(scale)
                v = v.multiply(scale)
                R1 *= scale
                R2 *= scale
            params = (p, v, R1, R2)

        if Stype == "plane":
            surfaces[Sid] = Plane(number, params, trsf)
        elif Stype == "sphere":
            surfaces[Sid] = Sphere(number, params, trsf)
        elif Stype == "cylinder" or Stype == "can":
            surfaces[Sid] = Cylinder(number, params, trsf, Stype == "can")
        elif Stype == "cylinder_elliptic" or Stype == "ecan":
            surfaces[Sid] = EllipticCylinder(number, params, trsf, Stype == "ecan")
        elif Stype == "cylinder_hyperbolic":
            surfaces[Sid] = HyperbolicCylinder(number, params, trsf)
        elif Stype == "cone" or Stype == "tcone":
            surfaces[Sid] = Cone(number, params, trsf, Stype == "tcone")
        elif Stype == "cone_elliptic":
            surfaces[Sid] = EllipticCone(number, params, trsf)
        elif Stype == "hyperboloid":
            surfaces[Sid] = Hyperboloid(number, params, trsf)
        elif Stype == "ellipsoid":
            surfaces[Sid] = Ellipsoid(number, params, trsf)
        elif Stype == "paraboloid":
            surfaces[Sid] = Paraboloid(number, params, trsf)
        elif Stype == "torus":
            surfaces[Sid] = Torus(number, params, trsf)
        elif Stype == "box":
            surfaces[Sid] = Box(number, params, trsf)
        else:
            print("Undefined", Sid, Stype)
            print(MCNPtype, number, MCNPparams)

    return surfaces


def pointsToCoeffs(scf):
    # mcnp implementation to convert 3 point plane to
    # plane parameters

    tpp = [0] * 4
    for i in range(1, 4):
        j = i % 3 + 1
        k = 6 - i - j
        k -= 1
        j -= 1
        tpp[i - 1] = (
            scf[j] * (scf[k + 3] - scf[k + 6])
            + scf[j + 3] * (scf[k + 6] - scf[k])
            + scf[j + 6] * (scf[k] - scf[k + 3])
        )
        tpp[3] += scf[i - 1] * (scf[j + 3] * scf[k + 6] - scf[j + 6] * scf[k + 3])

    xm = 0
    coeff = [0] * 4
    for i in range(1, 5):
        if xm == 0 and tpp[4 - i] != 0:
            xm = 1 / tpp[4 - i]
        coeff[4 - i] = tpp[4 - i] * xm

    # coeff [0:3] a,b,c plane parameters
    # coeff [3]   d plane parameter
    # normalization is d set to one if origin is not in the plane


def get_parabola_parameters(eVal, eVect, T, U):
    iaxis, comp = U[1]
    center = FreeCAD.Vector(T)
    axis = FreeCAD.Vector(eVect[iaxis][0])
    e1 = eVal[(iaxis + 1) % 3]
    focal = comp / (4 * e1)
    if focal < 0:
        focal = -focal
        axis = -axis
    return (center, axis, focal)


def get_cylinder_parameters(eVal, eVect, T, k, iaxis):

    other1 = (iaxis + 1) % 3
    other2 = (iaxis + 2) % 3

    eMin = eVal[other1]
    eMaj = eVal[other2]

    axis = FreeCAD.Vector(np.transpose(eVect)[iaxis])
    pos = FreeCAD.Vector(T)
    if abs(eMin - eMaj) < 1.0e-5:
        radius = float(np.sqrt(k / eMaj))
        return "cylinder", (pos, axis, radius)
    else:
        iMin = other1
        iMaj = other2
        if abs(eMin) < abs(eMaj):
            eMin, eMaj = eMaj, eMin
            iMin, iMaj = iMaj, iMin

        majorRad = float(np.sqrt(abs(k / eMaj)))
        minorRad = float(np.sqrt(abs(k / eMin)))
        minorAxis = FreeCAD.Vector(eVect.T[iMin])  # define axis in global geometry
        majorAxis = FreeCAD.Vector(eVect.T[iMaj])
        if np.sign(eMaj) == np.sign(eMin):
            return "cylinder_elliptic", (
                pos,
                axis,
                [minorRad, majorRad],
                [minorAxis, majorAxis],
            )
        else:
            if np.sign(k) == np.sign(eMaj):
                return "cylinder_hyperbolic", (
                    pos,
                    axis,
                    [minorRad, majorRad],
                    [minorAxis, majorAxis],
                )
            else:
                return "cylinder_hyperbolic", (
                    pos,
                    axis,
                    [majorRad, minorRad],
                    [majorAxis, minorAxis],
                )


def get_cone_parameters(eVal, eVect, T, iaxis):

    other1 = (iaxis + 1) % 3
    other2 = (iaxis + 2) % 3
    pos = FreeCAD.Vector(T)

    if abs(eVal[other1] - eVal[other2]) < 1e-5:
        axis = FreeCAD.Vector(np.transpose(eVect)[iaxis])
        tan = float(np.sqrt(-eVal[other1] / eVal[iaxis]))
        return "cone", (pos, axis, tan, True)
    else:
        for i in range(3):
            if np.sign(eVal[(i + 1) % 3]) == np.sign(eVal[(i + 2) % 3]):
                iaxis = i
                other1 = (iaxis + 1) % 3
                other2 = (iaxis + 2) % 3
                break

        axis = FreeCAD.Vector(np.transpose(eVect)[iaxis])
        minAxis = FreeCAD.Vector(np.transpose(eVect)[other1])
        majAxis = FreeCAD.Vector(np.transpose(eVect)[other2])
        Ra = abs(1 / eVal[iaxis])
        Rmin = abs(1 / eVal[other1])
        Rmaj = abs(1 / eVal[other2])

        if Rmin > Rmaj:
            Rmin, Rmaj = Rmaj, Rmin
            minAxis, majAxis = majAxis, minAxis

        return "cone_elliptic", (pos, axis, Ra, [Rmin, Rmaj], [minAxis, majAxis], True)


def get_hyperboloid_parameters(eVal, eVect, T, k, iaxis):
    cylTan = 1e3
    coneRad = 0.1

    ellipsoid = False
    if iaxis is None:
        iaxis = np.argmin(np.abs(eVal))
        ellipsoid = True

    other1 = (iaxis + 1) % 3
    other2 = (iaxis + 2) % 3
    Rad1 = float(np.sqrt(abs(k / eVal[other1])))
    Rad2 = float(np.sqrt(abs(k / eVal[other2])))
    other = other1 if Rad1 > Rad2 else other2

    majorRad = float(np.sqrt(abs(k / eVal[iaxis])))
    minorRad = float(np.sqrt(abs(k / eVal[other])))
    oneSheet = np.sign(k) != np.sign(eVal[iaxis])

    axis = FreeCAD.Vector(np.transpose(eVect)[iaxis])
    pos = FreeCAD.Vector(T)
    minorAxis = FreeCAD.Vector(eVect.T[other])  # define axis in global geometry
    majorAxis = FreeCAD.Vector(eVect.T[iaxis])

    t = majorRad / minorRad

    if t > cylTan and oneSheet:
        return get_cylinder_parameters(eVal, eVect, T, k, iaxis)
    elif minorRad < coneRad:
        return get_cone_parameters(eVal, eVect, T, iaxis)
    else:
        if ellipsoid:
            print("ellipical hyperboloid not implemented")
            print("single radius from {} eigen Value will be used".format(minorRad))
        return "hyperboloid", (
            pos,
            axis,
            [minorRad, majorRad],
            [minorAxis, majorAxis],
            oneSheet,
        )


def get_ellipsoid_parameters(eVal, eVect, T, k):

    cylTan = 1e3
    iaxis = None
    for i in range(3):
        if abs(eVal[i] - eVal[(i + 1) % 3]) < 1e-3:
            iaxis = (i + 2) % 3
            break

    if iaxis is None:
        print(eVal)
        print("cannot produce three radii Ellipsoid")
        return "ellipsoid_general", None

    other1 = (iaxis + 1) % 3
    iMaj = iaxis
    iMin = other1

    eMaj = eVal[iaxis]
    eMin = eVal[other1]

    if eMin < eMaj:
        eMin, eMaj = eMaj, eMin
        iMin, iMaj = iMaj, iMin

    RMaj = float(np.sqrt(abs(k / eMaj)))
    RMin = float(np.sqrt(abs(k / eMin)))
    majorAxis = FreeCAD.Vector(np.transpose(eVect)[iMaj])
    minorAxis = FreeCAD.Vector(np.transpose(eVect)[iMin])
    pos = FreeCAD.Vector(T)

    t = RMaj / RMin
    if t > cylTan:
        return get_cylinder_parameters(eVal, eVect, T, k, iMaj)
    else:
        return "ellipsoid", (pos, iaxis, [RMin, RMaj], [minorAxis, majorAxis])


def getGQAxis(eVal, k):

    # check if there is two equal eigenValues
    iaxis = None
    for i in range(3):
        if abs(eVal[i] - eVal[(i + 1) % 3]) < 1e-5:
            iaxis = (i + 2) % 3
            break

    if iaxis is None:
        iaxis = np.argmin(np.abs(eVal))

    e0 = eVal[iaxis]
    e1 = eVal[(iaxis + 1) % 3]
    e2 = eVal[(iaxis + 2) % 3]

    if k == 0:  # k == 0
        if (
            e0 == 0
        ):  # e1*X^2 + e2*Y^2             = 0    Intersecting  planes (real or imaginary)
            ek = None
        elif np.sign(e0) == np.sign(e1) and np.sign(e1) == np.sign(
            e2
        ):  # e1*X^2 + e2*Y^2 + e0*Z^2    = 0    Imaginary ellipsoid
            ek = None
        else:  # e1*X^2 + e2*Y^2 - e0*Z^2    = 0    Elliptic Cone
            ek = (-1, 0)

    elif np.sign(k) == np.sign(e1):  # e1 and k same sign  (e1 > 0)
        if e0 == 0:
            if np.sign(e1) == np.sign(
                e2
            ):  # e1*X^2 + e2*Y^2          + |k| = 0  Imaginary Elliptic cylinder
                ek = None
            else:  # e1*X^2 - e2*Y^2          + |k| = 0  Hyperpolic cylinder
                ek = (0, -1)
        elif np.sign(e0) == np.sign(e1) and np.sign(e1) == np.sign(
            e2
        ):  # e1*X^2 + e2*Y^2 + e0*Z^2 + |k| = 0  Imaginary ellipsoid
            ek = None
        else:  # e1*X^2 + e2*Y^2 - e0*Z^2 + |k| = 0  Hyperboloid
            ek = (-1, 1)

    else:  # e1 and k different sign
        if e0 == 0:  # e1*X^2 + e2*Y^2          - |k| = 0  Elliptic cylinder
            ek = (0, -1)
        elif np.sign(e0) == np.sign(e1) and np.sign(e1) == np.sign(
            e2
        ):  # e1*X^2 + e2*Y^2 + e0*Z^2 - |k| = 0  Ellipsoid
            ek = (1, -1)
        else:  # e1*X^2 + e2*Y^2 - e0*Z^2 - |k| = 0  Hyperbpoloid
            ek = (-1, 1)

    return iaxis, ek


def sq2params(params):

    a, b, c, d, e, f, g = params[0:7]
    x, y, z = params[7:10]

    gqa, gqb, gqc = a, b, c
    gqd, gqe, gqf = 0, 0, 0
    gqg = 2 * (d - a * x)
    gqh = 2 * (e - b * y)
    gqj = 2 * (f - c * z)
    gqk = g + a * x * x + b * y * y + c * z * z - 2 * d * x - 2 * e * y - 2 * f * z

    return gq2params((gqa, gqb, gqc, gqd, gqe, gqf, gqg, gqh, gqj, gqk))


def gq2params(x):
    # intial matrix: A
    #      a f g u
    #      f b h v
    #      g h c w
    #      u v w d

    # matrix displacement : D
    #      1 0 0 alpha
    #      0 1 0 beta
    #      0 0 1 gamma
    #      0 0 0 1

    #   Transpose(D) * A * D = mat3
    #   matrix mat3
    #      a f g 0
    #      f b h 0
    #      g h c 0
    #      0 0 0 (d + u*alpha + v*beta + w*gamma)

    #  with vector(alpha,beta,gamma) solution of
    #      a f h     alpha      -u
    #      f b g  *  beta   =   -v
    #      h g c     gamma      -w

    mat3 = np.array(
        [
            [x[0], x[3] / 2, x[5] / 2],
            [x[3] / 2, x[1], x[4] / 2],
            [x[5] / 2, x[4] / 2, x[2]],
        ]
    )
    X = np.array((x[6] / 2, x[7] / 2, x[8] / 2))

    # mat4 = np.array( [[x[0],x[3]/2,x[5]/2,x[6]/2], \
    #                  [x[3]/2,x[1],x[4]/2,x[7]/2], \
    #                  [x[5]/2,x[4]/2,x[2],x[8]/2], \
    #                  [x[6]/2,x[7]/2,x[8]/2,x[9]]] )

    # eigenValues and Vector
    eVal, vect = LA.eigh(mat3)
    XD = np.matmul(X, vect)  # X in diagonalised base

    Dinv = np.where(
        abs(eVal) < 1e-8, eVal, 1 / eVal
    )  # get inverse eigen value where eigen< 1e-8
    zero = (
        abs(eVal) < 1e-8
    ).nonzero()  # index in eigen value vector where eigen < 1e-8
    TD = -XD * Dinv  # Translation vector in diagonalized base

    k = np.matmul(TD, XD) + x[9]

    if len(zero) != 0:
        iz = zero[0]
        comp = 2 * XD[iz]
        if (
            abs(comp) > 1e-6
        ):  # zero eigenvalue but corresponding component in XD vector is non zero => paraboloid Curve => the k/comp value is the translation in this component direction
            TD[iz] = -k / comp
            U = (k, (iz, comp))
        else:
            U = (k, None)
    else:
        U = (k, None)

    T = np.matmul(TD, vect.T)

    return conicSurface(eVal, vect, T, U)


def conicSurface(eVal, vect, T, U):

    # paraboloid
    if U[1] is not None:
        params = get_parabola_parameters(eVal, vect, T, U)
        stype = "paraboloid"
        return stype, params

    # other conics
    k = U[0]
    iaxis, ek = getGQAxis(eVal, k)
    if ek == (0, -1):
        # Cylinder
        stype, params = get_cylinder_parameters(eVal, vect, T, -k, iaxis)
    elif ek == (-1, 0):
        # cone
        stype, params = get_cone_parameters(eVal, vect, T, iaxis)
    elif ek == (-1, 1):
        # hyperboloid
        # function can return cylinder or cone if hyperboiloid can be aproximated by such surfaces (within some criterion)
        stype, params = get_hyperboloid_parameters(eVal, vect, T, -k, iaxis)
    elif ek == (1, -1):
        # ellipsoid
        # function can return cylinder or cone if hyperboiloid can be aproximated by such surfaces (within some criterion)
        stype, params = get_ellipsoid_parameters(eVal, vect, T, -k)
    else:
        stype = "unknown"
        params = None
        print("No 2nd order or real surfaces")
    return stype, params
