import math
import os
import re
import xml.etree.ElementTree as ET

import FreeCAD
import numpy as np
from numpy import linalg as LA

from .Objects import CadCell, Cone, Cylinder, Plane, Sphere, Torus
from .XMLParser import get_cards


class XmlInput:
    def __init__(self, name):
        if not os.path.isfile(name):
            raise FileNotFoundError(f"File {name} does not exist")

        tree = ET.parse(name)
        root = tree.getroot()

        self.__inputcards__ = list(get_cards(root))
        return

    def GetFilteredCells(self, Ustart, depth, matcel_list, settings):

        if depth == 0:
            Ukeys = (Ustart,)
        else:
            subUniverses = getSubUniverses(Ustart, self.Universes)
            subUniverses.add(Ustart)

            for lev, Univ in self.levels.items():
                if Ustart in Univ:
                    break
            else:
                raise ValueError(f"Universe {Ustart} not found in the model")

            if depth == -1:
                levelMax = len(self.levels) - 1
            else:
                levelMax = min(lev + depth, len(self.levels) - 1)

            levelUniverse = set()
            for clev in range(lev, levelMax + 1):
                for U in self.levels[clev]:
                    levelUniverse.add(U)

            # select only universes wich are subuniverse of Ustart
            subUniverses = subUniverses.intersection(levelUniverse)
            Ukeys = list(self.Universes.keys())
            for U in Ukeys:
                if U not in subUniverses:
                    Ukeys.remove(U)

        FilteredCells = {}
        for U in Ukeys:
            FilteredCells[U] = selectCells(self.Universes[U], matcel_list)
            processSurfaces(FilteredCells[U], self.surfaces)

        # change the surface name in surface dict
        newSurfaces = {}
        for k in self.surfaces.keys():
            newkey = self.surfaces[k].id
            newSurfaces[newkey] = self.surfaces[k]

        for U, universe in FilteredCells.items():
            # set cell as CAD cell Object
            for cname, c in universe.items():
                universe[cname] = CadCell(c, settings=settings)

        return FilteredCells, newSurfaces

    def GetLevelStructure(self):
        containers = []
        Universe_dict = {}
        containers_label = set()

        for c in self.__inputcards__:
            if c.type != "cell":
                continue
            if c.U not in Universe_dict.keys():
                Universe_dict[c.U] = {}
            Universe_dict[c.U].update({c.name: c})

            if c.FILL:
                containers.append(c)
                containers_label.append(c.FILL)

        if 0 in Universe_dict.keys():
            root_universe = 0
        else:
            for k in Universe_dict.keys():
                if k not in containers_label:
                    root_universe = k
                    break

        # check all Universe have container cell
        for k in Universe_dict.keys():
            if k not in containers_label and k != root_universe:
                raise RuntimeError(f"Universe {k} has not container cell.")

        currentLevel = [root_universe]
        nextLevel = []
        univLevel = {0: {root_universe}}
        level = 0

        while True:
            level += 1
            univLevel[level] = set()
            for c in reversed(containers):
                if c.U in currentLevel:
                    c.Level = level
                    nextLevel.append(c.FILL)
                    univLevel[level].add(c.FILL)
                    containers.remove(c)

            if nextLevel == []:
                break
            currentLevel = nextLevel
            nextLevel = []

        lmax = len(univLevel)
        del univLevel[lmax - 1]
        for k in univLevel.keys():
            univLevel[k] = tuple(univLevel[k])
        self.levels = univLevel
        self.Universes = Universe_dict
        return

    def GetCell(self, name, settings):
        for c in self.__inputcards__:
            if c.ctype != "cell":
                continue
            if c.name != name:
                continue

            processSurfaces({c.name: c}, self.surfaces)
            newSurfaces = {}
            for k in self.surfaces.keys():
                newkey = self.surfaces[k].id
                newSurfaces[newkey] = self.surfaces[k]

            c = CadCell(c, settings=settings)
            c.setSurfaces(newSurfaces)
            return c

    def GetCells(self, U=None, Fill=None):
        cell_cards = {}
        for c in self.__inputcards__:
            if c.type != "cell":
                continue
            U_cell = c.U
            Fill_cell = c.FILL
            if U is None and Fill is None:
                cell_cards[c.name] = c
            elif U_cell == U and U is not None:
                cell_cards[c.name] = c
            elif Fill_cell == Fill and Fill is not None:
                cell_cards[c.name] = c

        return cell_cards

    def GetSurfaces(self):
        surf_cards = {}
        number = 1
        scale = 1.0  # don't change CAD units
        for c in self.__inputcards__:
            if c.type != "surface":
                continue
            surf_cards[c.name] = (c.stype, c.scoefs, number)
            number += 1

        self.surfaces = Get_primitive_surfaces(surf_cards, scale)


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
                    selected[name] = c  # Fill cell are not tested against material number
        elif config["cell"][0] == "exclude":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT not in config["mat"][1]:
                        if name not in config["cell"][1]:
                            selected[name] = c
                else:
                    if name not in config["cell"][1]:
                        selected[name] = c  # Fill cell are not tested against material number
        elif config["cell"][0] == "include":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT not in config["mat"][1]:
                        if name in config["cell"][1]:
                            selected[name] = c
                else:
                    if name in config["cell"][1]:
                        selected[name] = c  # Fill cell are not tested against material number

    # options are 'include' material
    elif config["mat"][0] == "include":
        if config["cell"][0] == "all":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT in config["mat"][1]:
                        selected[name] = c
                else:
                    selected[name] = c  # Fill cell are not tested against material number
        elif config["cell"][0] == "exclude":
            for c in cellList:
                if c.FILL is None:
                    if c.MAT in config["mat"][1]:
                        if name not in config["cell"][1]:
                            selected[name] = c
                else:
                    if name not in config["cell"][1]:
                        selected[name] = c  # Fill cell are not tested against material number
        elif config["cell"][0] == "include":
            for name, c in cellList.items():
                if c.FILL is None:
                    if c.MAT in config["mat"][1]:
                        if name in config["cell"][1]:
                            selected[name] = c
                else:
                    if name in config["cell"][1]:
                        selected[name] = c  # Fill cell are not tested against material number

    # remove complementary in cell of the universe
    # for cname,c in selected.items() :
    #   c.geom = remove_hash(cellList,cname)

    if not selected:
        raise ValueError("No cells selected. Check input or selection criteria in config file.")

    return selected


def processSurfaces(UCells, Surfaces):
    number = re.compile(r"\#?\s*\d+")

    for cname, c in UCells.items():
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

    surfaces = {}
    for Sid in mcnp_surfaces.keys():
        MCNPtype = mcnp_surfaces[Sid][0]
        MCNPparams = mcnp_surfaces[Sid][1]
        number = mcnp_surfaces[Sid][2]

        params = []
        Stype = None
        if MCNPtype in ("plane", "x-plane", "y-plane", "z-plane"):
            Stype = "plane"
            if MCNPtype == "plane":
                normal = FreeCAD.Vector(MCNPparams[0:3])
                params = (normal, MCNPparams[3] * scale)
            elif MCNPtype == "x-plane":
                params = (X_vec, MCNPparams[0] * scale)
            elif MCNPtype == "y-plane":
                params = (Y_vec, MCNPparams[0] * scale)
            elif MCNPtype == "z-plane":
                params = (Z_vec, MCNPparams[0] * scale)

        elif MCNPtype == "sphere":
            Stype = "sphere"
            params = (FreeCAD.Vector(MCNPparams[0:3]) * scale, MCNPparams[3] * scale)

        elif MCNPtype in ("x-cylinder", "y-cylinder", "z-cylinder"):
            R = MCNPparams[2]
            x1 = MCNPparams[0]
            x2 = MCNPparams[1]
            Stype = "cylinder"

            if MCNPtype == "x-cylinder":
                v = X_vec
                p = FreeCAD.Vector(0.0, x1, x2)
            elif MCNPtype == "y-cylinder":
                v = Y_vec
                p = FreeCAD.Vector(x1, 0.0, x2)
            elif MCNPtype == "z-cylinder":
                v = Z_vec
                p = FreeCAD.Vector(x1, x2, 0.0)

            if scale != 1.0:
                p = p.multiply(scale)
                R *= scale

            params = (p, v, R)

        elif MCNPtype in ("x-cone", "y-cone", "z-cone"):
            Stype = "cone"
            x1 = MCNPparams[0]
            x2 = MCNPparams[1]
            x3 = MCNPparams[2]
            p = FreeCAD.Vector(x1, x2, x3)
            t2 = MCNPparams[3]
            t = math.sqrt(t2)
            dblsht = True

            if MCNPtype == "x-cone":
                v = X_vec
            elif MCNPtype == "y-cone":
                v = Y_vec
            elif MCNPtype == "z-cone":
                v = Z_vec

            p = p.multiply(scale)
            params = (p, v, t, dblsht)

        elif MCNPtype in ["x-torus", "y-torus", "z-torus"]:
            Stype = "torus"
            p = FreeCAD.Vector(MCNPparams[0:3])
            Ra, r1, r2 = MCNPparams[3:6]

            if MCNPtype == "x-torus":
                v = X_vec
            elif MCNPtype == "y-torus":
                v = Y_vec
            elif MCNPtype == "z-torus":
                v = Z_vec

            if scale != 1.0:
                Ra *= scale
                r1 *= scale
                r2 *= scale
                p = p.multiply(scale)

            params = (p, v, Ra, r1, r2)

        elif MCNPtype == "quadric":
            Qparams = tuple(MCNPparams[0:10])
            Stype, quadric = gq2cyl(Qparams)

            if Stype == "cylinder":
                p = FreeCAD.Vector(quadric[0:3])
                v = FreeCAD.Vector(quadric[3:6])
                R = quadric[6]
                if scale != 1.0:
                    R *= scale
                    p = p.multiply(scale)

                params = (p, v, R)

            elif Stype == "cone":
                p = FreeCAD.Vector(quadric[0:3])
                v = FreeCAD.Vector(quadric[3:6])
                t = quadric[6]
                dblsht = quadric[7]
                if scale != 1.0:
                    p = p.multiply(scale)

                params = (p, v, t, dblsht)

            else:
                print(Stype)
                params = None
        #                get_quadric_surface(params)

        if Stype == "plane":
            surfaces[Sid] = Plane(Sid, number, params)
        elif Stype == "sphere":
            surfaces[Sid] = Sphere(Sid, number, params)
        elif Stype == "cylinder":
            surfaces[Sid] = Cylinder(Sid, number, params)
        elif Stype == "cone":
            surfaces[Sid] = Cone(Sid, number, params)
        elif Stype == "torus":
            surfaces[Sid] = Torus(Sid, number, params)
        else:
            print("Undefined", Sid)
            print(MCNPtype, number, MCNPparams)

    return surfaces


def gq2cyl(x):
    # Conversion de GQ a Cyl
    # Ax2+By2+Cz2+Dxy+Eyz+Fxz+Gx+Hy+Jz+K=0
    # x.T*M*x + b.T*x + K = 0
    minWTol = 5.0e-2
    minRTol = 1.0e-3
    # minRTol=3.e-1
    # lx = np.array(x)
    tp = ""
    M = np.array(
        [
            [x[0], x[3] / 2, x[5] / 2],
            [x[3] / 2, x[1], x[4] / 2],
            [x[5] / 2, x[4] / 2, x[2]],
        ]
    )
    w, P = LA.eigh(M)
    sw = np.sort(w)
    aw = np.abs(w)
    asw = np.sort(aw)
    # Test for cylinder (least abs value is much less than others)
    if asw[0] < minWTol * asw[1]:
        tp = "cylinder"
        rv = [0.0] * 7  # X0,Y0,Z0, VX, VY, VZ, R
        iaxis = np.where(aw == asw[0])[0][0]
        otherAxes = ((iaxis + 1) % 3, (iaxis + 2) % 3)
        if abs(w[otherAxes[0]] - w[otherAxes[1]]) > minRTol * asw[2]:
            tp = "not found - ellipsoid cylinder"
            rv = [0]
            return tp, rv
        # Vector de desplazamiento
        # x0 = -0.5*Pt*D-1*P*b pero ojo que un lambda es cero
        # P es la matriz de valores propios
        b = np.array(x[6:9])
        Pb = np.matmul(P, b)
        for i in otherAxes:
            Pb[i] /= w[i]
        x0 = -0.5 * np.matmul(P.T, Pb)
        k = -0.5 * np.matmul(x0, b) - x[9]
        # Resultados finales

        rv[0:3] = x0  # Punto del eje
        rv[3:6] = P[:, iaxis]  # Vector director
        rv[6] = np.sqrt(k / sw[1])  # Radio
    # Test for cone (incomplete, returns empty data list)
    elif np.sign(sw[0]) != np.sign(sw[2]):  # maybe cone
        tp = "cone"
        rv = [0.0] * 8  #  X0, Y0, Z0, VX, VY, VZ, tgAlpha, double sheet
        if np.sign(sw[0]) == np.sign(sw[1]):
            iaxis = np.where(w == sw[2])[0][0]
        else:
            iaxis = np.where(w == sw[0])[0][0]
        otherAxes = ((iaxis + 1) % 3, (iaxis + 2) % 3)
        if abs(w[otherAxes[0]] - w[otherAxes[1]]) > minRTol * asw[2]:
            tp = "not found - ellipsoid cone/hyperboloid"
            rv = [0]
            return tp, rv
        # Displacement vector ( x0 = -0.5*M^-1*b = -0.5*P.T*D^-1*P*b
        b = np.array(x[6:9])
        x0 = -0.5 * np.matmul(P, np.matmul(P.T, b) / w)
        k = x0.T @ M @ x0 - x[9]
        if np.abs(k * w[iaxis]) > minRTol * minRTol * asw[2]:

            # tp = 'not found - hyperboloid'
            # rv = [0]
            # force cone surface
            print("Force cone surface")
            tp = "cone"
            rv[0:3] = x0  # vertex point
            rv[3:6] = P[:, iaxis]  # axis direction
            rv[6] = np.sqrt(-w[otherAxes[0]] / w[iaxis])  # semiangle tangent
            rv[7] = True  # here always double sheet cones
            return tp, rv
        # return value
        rv[0:3] = x0  # vertex point
        rv[3:6] = P[:, iaxis]  # axis direction
        rv[6] = np.sqrt(-w[otherAxes[0]] / w[iaxis])  # semiangle tangent
        rv[7] = True  # here always double sheet cones
    else:
        tp = "not found - unknown"
        rv = [0]
    return tp, rv
