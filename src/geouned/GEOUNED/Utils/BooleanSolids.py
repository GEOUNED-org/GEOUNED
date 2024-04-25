#
#   Conversion to MCNP v0.0
#   Only one solid and planar surfaces
#
import math

import FreeCAD

from ..Utils.booleanFunction import BoolSequence
from ..Utils.Functions import GEOUNED_Surface, splitBOP
from ..Utils.Options.Classes import Options as opt

BoolVals = (None, True, False)


class CTelement:
    def __init__(self, val=None, S1=None, S2=None):
        self.diagonal = False
        self.S1 = S1
        self.S2 = S2
        if val is None:
            self.val = None
        else:
            if type(val) is int:
                self.diagonal = True
                self.type = None
            else:
                self.type = val.count(0)
            self.val = val

    def getTranspose(self):
        if self.diagonal:
            return self.val
        return CTelement(
            (self.val[0], self.val[3], self.val[2], self.val[1]), self.S2, self.S1
        )

    def getDependence(self):
        if self.diagonal:
            if self.val == 0:
                return True, False
            elif self.val == -1:
                return "Null", False
            else:
                return True, "Null"

        Ones = sum(self.val)
        if Ones == 4:
            return None, None
        elif Ones == 3:
            ind = self.val.index(0)
            if ind == 0:
                return False, None
            elif ind == 1:
                return True, None
            elif ind == 2:
                return None, True
            else:
                return None, False
        elif Ones == 2:
            ind1 = self.val.index(0)
            ind2 = self.val.index(0, ind1 + 1)
            if ind1 == 0 and ind2 == 1:
                return "Null", None
            elif ind1 == 2 and ind2 == 3:
                return None, "Null"
            elif ind1 == 0 and ind2 == 3:
                return False, False
            else:
                return True, True
        elif Ones == 1:
            ind = self.val.index(1)
            if ind == 0:
                return True, "Null"
            elif ind == 1:
                return False, "Null"
            elif ind == 2:
                return "Null", False
            else:
                return "Null", True
        else:
            return "Null", "Null"


class ConstraintTable(dict):

    def __init__(self):
        self.diagonal = None

    def __str__(self):

        varName = list(self.keys())

        # Constraint Table only diagonal terms
        if len(self[varName[0]]) == 1:
            line = ""
            for name in varName:
                element = self[name][name]
                line += " {:4d} : {}\n".format(name, element.val)
            return line

        outstr = "  "
        for name in varName:
            outstr = outstr + " {:3d}".format(name)
        outstr = outstr + "\n"

        for name1 in varName:
            line = " {:3d} ".format(name1)
            linenot = "~{:3d} ".format(name1)
            for name2 in varName:
                elmt = self[name1][name2]
                if elmt.diagonal:
                    line += " {:>2d} ".format(elmt.val)
                    linenot += "    "
                else:
                    line += " {}{} ".format(elmt.val[0], elmt.val[1])
                    linenot += " {}{} ".format(elmt.val[3], elmt.val[2])
            outstr += line + "\n"
            outstr += linenot + "\n"
        return outstr

    def addElement(self, k1, k2, val):

        if k1 in self.keys():
            self[k1][k2] = val
        else:
            self[k1] = {k2: val}

    def fillMissingElements(self):
        keys = list(self.keys())
        missing = []
        for i, k1 in enumerate(keys):
            for k2 in keys[i + 1 :]:
                if k2 in self[k1].keys():
                    elmt = self[k1][k2]
                    self[k2][k1] = elmt.getTranspose()
                else:
                    missing.append((k1, k2))

        for k1, k2 in missing:
            diag1 = self[k1][k1]
            diag2 = self[k2][k2]
            new = combineDiagElements(diag1, diag2)
            new.S1 = k1
            new.S2 = k2
            self[k1][k2] = new
            self[k2][k1] = new.getTranspose()

    def getOutSurfaces(self):
        out = []
        for k in self.keys():
            if self[k][k].val != 0:
                out.append(k)
        return out

    def getConstraintSet(self, valname):
        trueSet = {}
        falseSet = {}
        TNull = False
        FNull = False
        for k in self.keys():
            TValue, FValue = self[valname][k].getDependence()
            if TValue == "Null":
                TNull = True
                TValue = None
            if FValue == "Null":
                FNull = True
                FValue = None
            if TValue is not None:
                trueSet[k] = TValue
            if FValue is not None:
                falseSet[k] = FValue

        if TNull:
            trueSet = None
        if FNull:
            falseSet = None
        return trueSet, falseSet

    def solidInBox(self, Seq):  #  Sequence of the cell
        surfs = tuple(Seq.getSurfacesNumbers())
        if self.diagonal:
            seqValues = dict()
            for s in surfs:
                val = BoolVals[self[s][s].val]
                if val is not None:
                    seqValues[s] = val
            # if evaluate   None  : Box intersection Cell != 0      Part of the cell in the box
            #               True  : Box intersection Cell == Box    Cell cover the full region of the box. Void cell doesn't exist
            #               False : Box intersection Cell == 0      Cell out of the box
            res = Seq.evaluate(seqValues)
            return res if type(res) is bool else None

        else:
            trueSet, falseSet = self.getConstraintSet(surfs[0])
            if trueSet is not None:
                trueVal = Seq.evaluate(trueSet)
                trueVal = trueVal if type(trueVal) is bool else None
                if trueVal is None:
                    return None
                if falseSet is not None:
                    falseVal = Seq.evaluate(falseSet)
                    falseVal = falseVal if type(falseVal) is bool else None
                    if falseVal is None:
                        return None  # Part of the cell in the box
                    if trueVal == falseVal:
                        return trueVal  # True Cover full cell, False not in the box
                    else:
                        return None  # Part of the cell in the box
                else:
                    return trueVal
            elif falseSet is not None:
                falseVal = Seq.evaluate(falseSet)
                return falseVal if type(falseVal) is bool else None
            else:
                print("Bad trouble surfaces {} is on none side of the box!!")
                return False


def combineDiagElements(d1, d2):

    if d1.val == 0 and d2.val == 0:
        return CTelement((1, 1, 1, 1))
    elif d1.val == 1 and d2.val == 0:
        return CTelement((1, 1, 0, 0))
    elif d1.val == -1 and d2.val == 0:
        return CTelement((0, 0, 1, 1))
    elif d1.val == 0 and d2.val == 1:
        return CTelement((1, 0, 0, 1))
    elif d1.val == 0 and d2.val == -1:
        return CTelement((0, 1, 1, 0))
    elif d1.val == 1 and d2.val == 1:
        return CTelement((1, 0, 0, 0))
    elif d1.val == 1 and d2.val == -1:
        return CTelement((0, 1, 0, 0))
    elif d1.val == -1 and d2.val == 1:
        return CTelement((0, 0, 0, 1))
    elif d1.val == -1 and d2.val == -1:
        return CTelement((0, 0, 1, 0))


def buildCTableFromSolids(Box, SurfInfo, option="diag"):

    if type(SurfInfo) is dict:
        surfaces = SurfInfo
        surfaceList = tuple(surfaces.keys())
    elif type(SurfInfo) is tuple:
        surfaceList, surfaces = SurfInfo
    else:
        surfaces = SurfInfo.Surfaces
        surfaceList = SurfInfo.surfaceList

    if type(surfaces[surfaceList[0]]) is GEOUNED_Surface:
        for s in surfaceList:
            ss = surfaces[s]
            ss.__boundBox__ = Box.BoundBox
            ss.buildSurface()
    else:
        for s in surfaceList:
            surfaces[s].buildShape(Box.BoundBox)

    CTable = ConstraintTable()
    if option == "diag":
        CTable.diagonal = True
    else:
        CTable.diagonal = False

    for i, s1 in enumerate(surfaceList):
        res, splitRegions = splitSolid_fast(Box, surfaces[s1], True)
        # res,splitRegions = splitSolid_fast(Box,Surfaces.getSurface(s1),True)

        CTable.addElement(s1, s1, CTelement(res, s1, s1))
        if option == "diag":
            continue
        if splitRegions is None:
            continue  # loop, no region to be split by s2
        for j, s2 in enumerate(surfaceList[i + 1 :]):
            posS1, negS1 = splitRegions

            pos0 = None
            for solid in posS1:

                pos = splitSolid_fast(solid, surfaces[s2], False)

                # pos = splitSolid_fast(solid,Surfaces.getSurface(s2),False)
                if pos == (1, 1):
                    break  # s2 intersect S1 Region
                if pos0 is None:
                    pos0 = pos
                else:
                    if pos != pos0:  # s1 regions are on both side of s2
                        pos = (1, 1)
                        break

            neg0 = None
            for solid in negS1:
                # neg = splitSolid_fast(solid,Surfaces.getSurface(s2),False)
                neg = splitSolid_fast(solid, surfaces[s2], False)
                if neg == (1, 1):
                    break  # s2 intersect S1 Region
                if neg0 is None:
                    neg0 = neg
                else:
                    if neg != neg0:  # s1 regions are on both side of s2
                        neg = (1, 1)
                        break

            val = (pos[0], pos[1], neg[1], neg[0])
            CTable.addElement(s1, s2, CTelement(val, s1, s2))

    # if some surfaces don't cross the box some elements in Constraint table are not filled
    if option != "diag":
        CTable.fillMissingElements()
    return CTable


def removeExtraSurfaces(CellSeq, CTable):
    # checking is make on solid cell definition to be removed from void cell
    outSurfaces = set(CTable.getOutSurfaces())
    newDef = BoolSequence(operator="OR")

    # Loop over all compound solids of the metaSolid

    if CellSeq.level == 0 and len(CellSeq.elements) == 1:
        return CellSeq

    if CellSeq.operator == "AND":
        newSeq = BoolSequence(operator="AND")
        newSeq.append(CellSeq.copy())
        CellSeq.assign(newSeq)

    for subCell in CellSeq.elements:
        nullcell = False
        subCell.check()
        if type(subCell.elements) is bool:
            chk = not subCell.elements
        else:
            chk = None

        if chk == False:  # the cell doesn't exist
            nullcell = True
        elif chk == True:  # the cell describe the full universe
            newDef.elements = True
            newDef.level = -1
            return newDef

        # if subcell has finite volume check it intersection with the box
        if not nullcell:
            res = CTable.solidInBox(subCell)
            if res == None:
                # subcell intersect the box
                # get the surfaces of the solids out of the box
                # get reduced definition

                # if subcell lev!= 0 remove surface operation is not valid
                if subCell.level == 0:
                    removeSurf = outSurfaces & subCell.getSurfacesNumbers()
                    for s in removeSurf:
                        val = True if CTable[s][s].val > 0 else False
                        subCell.substitute(s, val)

                if type(subCell.elements) is bool:
                    if subCell.elements == False:  #  cell does not intersect void box
                        continue
                    else:  # cell cover fully void box
                        newDef.elements = True
                        newDef.level = -1
                        return newDef
                else:
                    newDef.append(subCell)

            elif res == True:
                # subcell cover the full box region Void cell doesn't exist
                newDef.elements = True
                newDef.level = -1
                return newDef

    newDef.clean()

    return newDef


def splitSolid_fast(solid, surf, box):

    if box:
        if surf.shape:
            comsolid = splitBOP(solid, [surf.shape], opt.splitTolerance)
        else:
            return checkSign(solid, surf), None

        if len(comsolid.Solids) <= 1:
            return checkSign(solid, surf), None
        # sgn = checkSign(solid,surf)   # if "box" and single object => the box is not split, surface s1 out of the box.
        # if sgn == 1 :                 # The sign is the side of surface s1 where the box is located
        #    # return ((1,0),(0,0)),None  # return the diagonal element of the Constraint Table for s1
        #    return ((1,0),(0,0)),None  # return the diagonal element of the Constraint Table for s1
        # else:
        #    return ((0,0),(0,1)),None

        else:
            posSol = []
            negSol = []
            for s in comsolid.Solids:
                sgn = checkSign(s, surf)
                if sgn == 1:
                    posSol.append(s)
                else:
                    negSol.append(s)
            return 0, (
                posSol,
                negSol,
            )  # return the diagonal element of the Constraint Table for s1, and solids to be split by s2
            # return ((1,0),(0,1)), (posSol,negSol) # return the diagonal element of the Constraint Table for s1, and solids to be split by s2

    else:
        # "not box" => return the position of the +/- region of s1 (the solid) with respect s2
        if surf.shape:
            dist = solid.distToShape(surf.shape)[0]
        else:
            dist = 1.0
        if dist > 1e-6:  # face doesn't intersect solid
            sgn = checkSign(solid, surf)
            if sgn == 1:
                return (
                    1,
                    0,
                )  # values of s2 and -s2   "0" means this region doesn't exist
            else:
                return (0, 1)
        else:
            return (1, 1)  # values of s2 and -s2


# find one point inside a solid (region)
def point_inside(solid):

    point = solid.CenterOfMass
    if solid.isInside(point, 0.0, False):
        return point

    cut_line = 32
    cut_box = 2

    v1 = solid.Vertexes[0].Point
    for vi in range(len(solid.Vertexes) - 1, 0, -1):
        v2 = solid.Vertexes[vi].Point
        dv = (v2 - v1) * 0.5

        n = 1
        while True:
            for i in range(n):
                point = v1 + dv * (1 + 0.5 * i)
                if solid.isInside(point, 0.0, False):
                    return point
            n = n * 2
            dv = dv * 0.5
            if n > cut_line:
                break

    BBox = solid.optimalBoundingBox(False)
    box = [BBox.XMin, BBox.XMax, BBox.YMin, BBox.YMax, BBox.ZMin, BBox.ZMax]

    boxes, centers = CutBox(box)
    n = 0

    while True:
        for p in centers:
            pp = FreeCAD.Vector(p[0], p[1], p[2])
            if solid.isInside(pp, 0.0, False):
                return pp

        subbox = []
        centers = []
        for b in boxes:
            btab, ctab = CutBox(b)
            subbox.extend(btab)
            centers.extend(ctab)
        boxes = subbox
        n = n + 1

        if n == cut_box:
            break

    return pointFromSurface(solid)


def pointFromSurface(solid):

    for face in solid.Faces:
        parameters = face.ParameterRange
        u = (parameters[0] + parameters[1]) / 2
        v = (parameters[2] + parameters[3]) / 2
        pface = face.valueAt(u, v)
        normal = face.normalAt(u, v)

        d = 10
        pp = pface + d * normal
        while d > 1e-8:
            if solid.isInside(pp, 0.0, False):
                return pp
            d *= 0.5
            pp = pface + d * normal

        normal = -normal
        d = 10
        pp = pface + d * normal
        while d > 1e-8:
            if solid.isInside(pp, 0.0, False):
                return pp
            d *= 0.5
            pp = pface + d * normal

    print("Solid not found in bounding Box (Volume : {})".format(solid.Volume))
    return None


# divide a box into 8 smaller boxes
def CutBox(Box):
    xmid = (Box[1] + Box[0]) * 0.5
    ymid = (Box[3] + Box[2]) * 0.5
    zmid = (Box[5] + Box[4]) * 0.5

    b1 = (Box[0], xmid, Box[2], ymid, Box[4], zmid)
    p1 = (0.5 * (Box[0] + xmid), 0.5 * (Box[2] + ymid), 0.5 * (Box[4] + zmid))

    b2 = (xmid, Box[1], Box[2], ymid, Box[4], zmid)
    p2 = (0.5 * (xmid + Box[1]), 0.5 * (Box[2] + ymid), 0.5 * (Box[4] + zmid))

    b3 = (Box[0], xmid, ymid, Box[3], Box[4], zmid)
    p3 = (0.5 * (Box[0] + xmid), 0.5 * (ymid + Box[3]), 0.5 * (Box[4] + zmid))

    b4 = (xmid, Box[1], ymid, Box[3], Box[4], zmid)
    p4 = (0.5 * (xmid + Box[1]), 0.5 * (ymid + Box[3]), 0.5 * (Box[4] + zmid))

    b5 = (Box[0], xmid, Box[2], ymid, zmid, Box[5])
    p5 = (0.5 * (Box[0] + xmid), 0.5 * (Box[2] + ymid), 0.5 * (zmid + Box[5]))

    b6 = (xmid, Box[1], Box[2], ymid, zmid, Box[5])
    p6 = (0.5 * (xmid + Box[1]), 0.5 * (Box[2] + ymid), 0.5 * (zmid + Box[5]))

    b7 = (Box[0], xmid, ymid, Box[3], zmid, Box[5])
    p7 = (0.5 * (Box[0] + xmid), 0.5 * (ymid + Box[3]), 0.5 * (zmid + Box[5]))

    b8 = (xmid, Box[1], ymid, Box[3], zmid, Box[5])
    p8 = (0.5 * (xmid + Box[1]), 0.5 * (ymid + Box[3]), 0.5 * (zmid + Box[5]))

    return (b1, b2, b3, b4, b5, b6, b7, b8), (p1, p2, p3, p4, p5, p6, p7, p8)


def checkSign(solid, surf):

    point = point_inside(solid)

    if surf.Type == "Plane":
        r = point - surf.Surf.Position
        if surf.Surf.Axis.dot(r) > 0:
            return 1
        else:
            return -1

    elif surf.Type == "Cylinder":
        r = point - surf.Surf.Center
        L2 = r.Length * r.Length
        z = surf.Surf.Axis.dot(r)
        z2 = z * z
        R2 = surf.Surf.Radius * surf.Surf.Radius
        if L2 - z2 > R2:
            return 1
        else:
            return -1

    elif surf.Type == "Sphere":
        r = point - surf.Surf.Center
        if r.Length > surf.Surf.Radius:
            return 1
        else:
            return -1

    elif surf.Type == "Cone":
        r = point - surf.Surf.Apex
        r.normalize()
        z = round(surf.Surf.Axis.dot(r), 15)
        alpha = math.acos(z)

        if alpha > surf.Surf.SemiAngle:
            return 1
        else:
            return -1

    elif surf.Type == "Torus":
        r = point - surf.Surf.Center
        h = r.dot(surf.Surf.Axis)
        rho = r - h * surf.Surf.Axis

        rp = math.sqrt((rho.Length - surf.Surf.MajorRadius) ** 2 + h**2)
        if rp > surf.Surf.MinorRadius:
            return 1
        else:
            return -1
