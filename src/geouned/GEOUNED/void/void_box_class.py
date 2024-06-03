"""File with the VoidBox class"""

import logging

import FreeCAD
import Part

from ..conversion import cell_definition as Conv
from ..decompose import decom_one as Decom
from ..utils.basic_functions_part1 import is_opposite
from ..utils.boolean_function import BoolSequence
from ..utils.boolean_solids import build_c_table_from_solids, remove_extra_surfaces
from ..utils.functions import GeounedSolid, GeounedSurface

logger = logging.getLogger("general_logger")


class VoidBox:
    def __init__(self, MetaSolids, Enclosure):

        self.Objects = []
        if "BoundBox" in str(Enclosure):
            self.BoundBox = Enclosure
            self.PieceEnclosure = None
        else:
            self.BoundBox = Enclosure.BoundBox
            self.PieceEnclosure = Enclosure

        for m in MetaSolids:
            if not m.BoundBox:
                continue
            if m.BoundBox.isValid():
                if self.BoundBox.intersect(m.BoundBox):
                    Obj = self.copy_meta(m)
                    self.remove_extra_comp(Obj, self.BoundBox)
                    self.Objects.append(Obj)
        return

    def split(self, minSize=200):

        dims = [self.BoundBox.XLength, self.BoundBox.YLength, self.BoundBox.ZLength]
        coord = ["X", "Y", "Z"]
        for i in range(2, -1, -1):
            if 0.5 * dims[i] < minSize:
                del dims[i]
                del coord[i]

        if len(dims) == 0:
            return None

        ip = dims.index(max(dims))
        #    print ('dims : {} {}'.format(coord[ip],0.5*dims[ip]))
        if coord[ip] == "X":
            pos = self.BoundBox.XMin + 0.5 * self.BoundBox.XLength
            X1Min = self.BoundBox.XMin
            X1Max = pos
            X2Min = pos
            X2Max = self.BoundBox.XMax
            Y1Min = self.BoundBox.YMin
            Y1Max = self.BoundBox.YMax
            Y2Min = Y1Min
            Y2Max = Y1Max
            Z1Min = self.BoundBox.ZMin
            Z1Max = self.BoundBox.ZMax
            Z2Min = Z1Min
            Z2Max = Z1Max
        elif coord[ip] == "Y":
            pos = self.BoundBox.YMin + 0.5 * self.BoundBox.YLength
            X1Min = self.BoundBox.XMin
            X1Max = self.BoundBox.XMax
            X2Min = X1Min
            X2Max = X1Max
            Y1Min = self.BoundBox.YMin
            Y1Max = pos
            Y2Min = pos
            Y2Max = self.BoundBox.YMax
            Z1Min = self.BoundBox.ZMin
            Z1Max = self.BoundBox.ZMax
            Z2Min = Z1Min
            Z2Max = Z1Max
        else:
            pos = self.BoundBox.ZMin + 0.5 * self.BoundBox.ZLength
            X1Min = self.BoundBox.XMin
            X1Max = self.BoundBox.XMax
            X2Min = X1Min
            X2Max = X1Max

            Y1Min = self.BoundBox.YMin
            Y1Max = self.BoundBox.YMax
            Y2Min = Y1Min
            Y2Max = Y1Max

            Z1Min = self.BoundBox.ZMin
            Z1Max = pos
            Z2Min = pos
            Z2Max = self.BoundBox.ZMax

        VMin1 = FreeCAD.Vector(X1Min, Y1Min, Z1Min)
        VMax1 = FreeCAD.Vector(X1Max, Y1Max, Z1Max)
        VMin2 = FreeCAD.Vector(X2Min, Y2Min, Z2Min)
        VMax2 = FreeCAD.Vector(X2Max, Y2Max, Z2Max)
        box1 = FreeCAD.BoundBox(VMin1, VMax1)
        box2 = FreeCAD.BoundBox(VMin2, VMax2)

        if self.PieceEnclosure is None:
            Space1 = VoidBox(self.Objects, box1)
            Space2 = VoidBox(self.Objects, box2)
            VoidBoxTuple = (Space1, Space2)
        else:
            Space1 = self.piece_enclosure_split(box1)
            Space2 = self.piece_enclosure_split(box2)
            VoidBoxTuple = (Space1, Space2)
            if Space1 is None:
                VoidBoxTuple = (Space2,)
            if Space2 is None:
                VoidBoxTuple = (Space1,)
            if Space1 is None and Space2 is None:
                VoidBoxTuple = ()

        return VoidBoxTuple

    def piece_enclosure_split(self, Box, Tolerance=1.0e-13):
        """This function creates a box-shaped solid with the new limits of given bounding box and
        it is intersected with the piece of nested enclosure to create the new void cell.
        If the limited region does not intersect with the piece, no void cell is created.
        """

        Cube = Part.makeBox(
            Box.XLength,
            Box.YLength,
            Box.ZLength,
            FreeCAD.Vector(Box.XMin, Box.YMin, Box.ZMin),
            FreeCAD.Vector(0, 0, 1),
        )
        dist = Cube.distToShape(self.PieceEnclosure)[0]
        try:
            if abs(dist / Box.DiagonalLength) > Tolerance:
                return None
        except ZeroDivisionError:
            return None
        ShapeObject = Cube.common(self.PieceEnclosure)
        try:
            reldif = (Cube.Volume - ShapeObject.Volume) / Cube.Volume
        except ZeroDivisionError:
            return None
        if abs(reldif) <= Tolerance:
            return VoidBox(self.Objects, Box)
        elif ShapeObject.Solids:
            Solid = ShapeObject.Solids[0]
            return VoidBox(self.Objects, Solid)
        else:
            return None

    def refine(self):
        Cube = Part.makeBox(
            self.BoundBox.XLength,
            self.BoundBox.YLength,
            self.BoundBox.ZLength,
            FreeCAD.Vector(self.BoundBox.XMin, self.BoundBox.YMin, self.BoundBox.ZMin),
            FreeCAD.Vector(0, 0, 1),
        )

        for m in self.Objects:
            self.remove_extra_comp(m, Cube, mode="dist")
        return

    def get_void_complementary(self, Surfaces, options, tolerances, numeric_format, simplify="no"):
        if self.PieceEnclosure is None:
            boxDef = BoolSequence(operator="AND")
            center = self.BoundBox.Center
            bBox = self.BoundBox
            for p in self.get_bound_planes():
                id, exist = Surfaces.add_plane(p, options, tolerances, numeric_format, False)
                if exist:
                    s = Surfaces.get_surface(id)
                    if is_opposite(p.Surf.Axis, s.Surf.Axis):
                        id = -id
                if is_opposite(p.Surf.Axis, p.Surf.Position - center):
                    boxDef.elements.append(id)
                else:
                    boxDef.elements.append(-id)
                enclosure = False
            d = options.enlargeBox

        else:
            UniverseBox = self.PieceEnclosure.BoundBox
            TempPieceEnclosure = GeounedSolid(None, self.PieceEnclosure)
            comsolid, err = Decom.SplitSolid(
                Part.makeCompound(TempPieceEnclosure.Solids),
                UniverseBox,
                options,
                tolerances,
                numeric_format,
            )
            Surfaces.extend(
                Decom.extract_surfaces(
                    comsolid,
                    "All",
                    UniverseBox,
                    options,
                    tolerances,
                    numeric_format,
                    MakeObj=True,
                ),
                options,
                tolerances,
                numeric_format,
            )
            TempPieceEnclosure.update_solids(comsolid.Solids)
            Conv.cellDef(
                TempPieceEnclosure,
                Surfaces,
                UniverseBox,
                options,
                tolerances,
                numeric_format,
            )

            boxDef = TempPieceEnclosure.Definition
            bBox = self.PieceEnclosure.BoundBox
            enclosure = True
            d = max(options.enlargeBox, 2)

        Box = Part.makeBox(
            bBox.XLength + 2 * d,
            bBox.YLength + 2 * d,
            bBox.ZLength + 2 * d,
            FreeCAD.Vector(bBox.XMin - d, bBox.YMin - d, bBox.ZMin - d),
            FreeCAD.Vector(0, 0, 1),
        )

        voidSolidDef = BoolSequence(operator="OR")

        cellIn = []
        for m in self.Objects:
            voidSolidDef.append(m.Definition)
            if m.IsEnclosure:
                continue
            cellIn.append(m.__id__)

        # voidSolidDef.join_operators()

        if not voidSolidDef.elements:
            return (
                boxDef,
                None,
            )  #  Cell to get complementary are null => Void is only box definition

        # join all basic solids into one big meta Object
        # CAD solid representation is not needed because
        # here we are working with surfaces and void box

        complementary = BoolSequence(operator="AND")
        complementary.append(boxDef)
        if simplify != "no":
            surfList = voidSolidDef.get_surfaces_numbers()

            if enclosure:
                surfList.update(boxDef.get_surfaces_numbers())
            else:
                for s in boxDef.elements:
                    val = s > 0
                    voidSolidDef.substitute(abs(s), val)
                voidSolidDef.clean()
            if type(voidSolidDef.elements) is bool:
                res = voidSolidDef.elements
            else:
                res = None

            if enclosure or res is None:
                surfaceDict = {}
                for i in surfList:
                    surfaceDict[i] = Surfaces.get_surface(i)
                CTable = build_c_table_from_solids(Box, surfaceDict, option=simplify, options=options)
            else:
                if res is True:
                    return None, None
                else:
                    return boxDef, None

            newTemp = BoolSequence(operator="OR")

            if voidSolidDef.level == 0:
                if len(voidSolidDef.elements) == 1:
                    voidSolidDef.operator = "AND"
                cellVoid = BoolSequence(operator="OR")
                cellVoid.append(voidSolidDef)
                voidSolidDef = cellVoid

            for solDef in voidSolidDef.elements:
                newSolid = remove_extra_surfaces(solDef, CTable)
                if type(newSolid.elements) is not bool:
                    newTemp.append(newSolid)
                elif newSolid.elements is True:
                    return None, None

            voidSolidDef = newTemp

        else:
            if voidSolidDef.level == 0:
                if len(voidSolidDef.elements) == 1:
                    voidSolidDef.operator = "AND"
                cellVoid = BoolSequence(operator="OR")
                cellVoid.append(voidSolidDef)
                voidSolidDef = cellVoid

        if voidSolidDef.elements is True:
            return None, None
        elif voidSolidDef.elements is False or voidSolidDef.elements == []:
            return boxDef, None

        if voidSolidDef.level == 0:
            compSeq = voidSolidDef.get_complementary()
        else:

            if voidSolidDef.level == 1 and voidSolidDef.operator == "AND":
                compSeq = BoolSequence(operator="OR")
            else:
                compSeq = BoolSequence(operator="AND")

            for comp in voidSolidDef.elements:
                if simplify == "no":
                    comp.check()
                    if type(comp.elements) is bool:
                        chk = comp.elements
                    else:
                        chk = None

                    # solid in cover full Void cell volume  => Void cell doesn't exist
                    if chk is True:
                        logger.warning("void Cell should not exist")
                        return None, None

                    # solid cell is not in void cell Void cell volume  => doesn't contribute to void definition
                    elif chk is False:
                        continue

                pmoc = comp.get_complementary()
                compSeq.append(pmoc)

        if simplify == "full":
            if enclosure:
                complementary.append(compSeq)
                complementary.simplify(CTable)
            else:
                compSeq.simplify(CTable)
                complementary.append(compSeq)
        else:
            compSeq.simplify(None)
            complementary.simplify(None)
            complementary.append(compSeq)

        complementary.clean()
        complementary.level_update()

        if type(complementary.elements) is bool:
            return None, None
        else:
            return complementary, cellIn

    def get_numbers(self):
        ns = 0
        nb = 0

        for m in self.Objects:
            ns += len(m.Surfaces)
            nb += len(m.Definition.elements) if m.Definition.level > 0 else 1

        return ns, nb

    def get_bound_planes(self):
        Xmid = 0.5 * (self.BoundBox.XMin + self.BoundBox.XMax)
        Ymid = 0.5 * (self.BoundBox.YMin + self.BoundBox.YMax)
        Zmid = 0.5 * (self.BoundBox.ZMin + self.BoundBox.ZMax)
        LX = self.BoundBox.ZMin + self.BoundBox.XLength
        LY = self.BoundBox.ZMin + self.BoundBox.YLength
        LZ = self.BoundBox.ZMin + self.BoundBox.ZLength
        PXMin = GeounedSurface(
            (
                "Plane",
                (
                    FreeCAD.Vector(self.BoundBox.XMin, Ymid, Zmid),
                    FreeCAD.Vector(1, 0, 0),
                    LY,
                    LZ,
                ),
            ),
            self.BoundBox,
        )
        PXMax = GeounedSurface(
            (
                "Plane",
                (
                    FreeCAD.Vector(self.BoundBox.XMax, Ymid, Zmid),
                    FreeCAD.Vector(-1, 0, 0),
                    LY,
                    LZ,
                ),
            ),
            self.BoundBox,
        )
        PYMin = GeounedSurface(
            (
                "Plane",
                (
                    FreeCAD.Vector(Xmid, self.BoundBox.YMin, Zmid),
                    FreeCAD.Vector(0, 1, 0),
                    LZ,
                    LX,
                ),
            ),
            self.BoundBox,
        )
        PYMax = GeounedSurface(
            (
                "Plane",
                (
                    FreeCAD.Vector(Xmid, self.BoundBox.YMax, Zmid),
                    FreeCAD.Vector(0, -1, 0),
                    LZ,
                    LX,
                ),
            ),
            self.BoundBox,
        )
        PZMin = GeounedSurface(
            (
                "Plane",
                (
                    FreeCAD.Vector(Xmid, Ymid, self.BoundBox.ZMin),
                    FreeCAD.Vector(0, 0, 1),
                    LX,
                    LY,
                ),
            ),
            self.BoundBox,
        )
        PZMax = GeounedSurface(
            (
                "Plane",
                (
                    FreeCAD.Vector(Xmid, Ymid, self.BoundBox.ZMax),
                    FreeCAD.Vector(0, 0, -1),
                    LX,
                    LY,
                ),
            ),
            self.BoundBox,
        )

        return (PXMin, PXMax, PYMin, PYMax, PZMin, PZMax)

    def remove_extra_comp(self, Obj, Box, mode="box"):
        reducedSol = []
        reducedDef = BoolSequence(operator="OR")
        if not Obj.Solids:
            return
        # Compare Solid BoundBox (here Box is BoundBox Object)
        if mode == "box":
            for i, sol in enumerate(Obj.Solids):
                if sol.BoundBox.isValid():
                    if Box.intersect(sol.BoundBox):
                        reducedSol.append(sol)
                        reducedDef.append(Obj.Definition.elements[i])

        # Compare solid using distToshape (here Box is a Solid Cube object)
        else:
            for i, sol in enumerate(Obj.Solids):
                dist = Box.distToShape(sol)[0]
                if dist == 0:
                    reducedSol.append(sol)
                    reducedDef.append(Obj.Definition.elements[i])

        if len(reducedSol) < len(Obj.Solids):
            Obj.update_solids(reducedSol)
            Obj.set_definition(reducedDef)
        return

    def copy_meta(self, m):
        solidsCopy = m.Solids[:]
        facesCopy = m.Faces[:]
        Meta = GeounedSolid(m.__id__, solidsCopy)
        Meta.set_definition(m.Definition.copy())
        Meta.set_faces(facesCopy)
        if m.IsEnclosure:
            Meta.IsEnclosure = True
            Meta.EnclosureID = m.EnclosureID
        return Meta
