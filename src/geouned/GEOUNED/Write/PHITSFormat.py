##############################
# Module to write PHITS input #
##############################

# Mainly based on the MCNPFormat.py
# Some are modified for PHITS
# 1. Chosen the graveyard surface as the boundary of the outer void region (-1)
# 2. Unified the auto genarated void regions when envelope*_*_ and/or enclosure*_*_ are not setted
# 3. Applied MATERIAL DEFINITION block when dummyMat = True
# 4. Applied VOLUME DEFINITION block when volCARD = True
# 5. Eliminated the only MCNP related parts
# 6. Added some comments to remind

import re
from datetime import datetime

import FreeCAD

from ..CodeVersion import *
from ..Utils.BasicFunctions_part1 import isOposite, pointsToCoeffs
from ..Utils.Functions import SurfacesDict
from ..Utils.Options.Classes import Options as opt
from ..Write.Functions import (
    CellString,
    PHITSSurface,
    changeSurfSign,
    writePHITSCellDef,
)


class PhitsInput:
    def __init__(self, Meta, Surfaces, setting):
        self.Title = setting["title"]
        self.VolSDEF = setting["volSDEF"]
        self.VolCARD = setting["volCARD"]
        self.U0CARD = setting["UCARD"]
        self.DummyMat = setting["dummyMat"]
        self.Matfile = setting["matFile"]
        self.voidMat = setting["voidMat"]
        self.startCell = setting["startCell"]
        self.Cells = Meta
        self.Options = {"Volume": self.VolCARD, "Universe": self.U0CARD}

        self.StepFile = setting["stepFile"]
        if isinstance(self.StepFile, (tuple, list)):
            self.StepFile = "; ".join(self.StepFile)

        if self.Title == "":
            self.Title = self.StepFile

        self.__getSurfaceTable__()
        self.__simplifyPlanes__(Surfaces)

        self.Surfaces = self.__sortedSurfaces__(Surfaces)
        self.Materials = set()

        return

    def writePHITS(self, filename):
        print(f"write PHITS file {filename}")
        self.inpfile = open(filename, "w", encoding="utf-8")
        self.__write_PHITS_header__(filename)

        cHeader = """\
$
$ ##########################################################
$                  CELL DEFINITION
$ ##########################################################
$
[CELL]\n"""

        self.inpfile.write(cHeader)
        self.__write_PHITS_cell_block__()
        self.inpfile.write(" \n")

        surfaceHeader = """\
$
$ ##########################################################
$                  SURFACE DEFINITION
$ ##########################################################
$
[SURFACE]\n"""

        self.inpfile.write(surfaceHeader)
        self.__write_PHITS_surface_block__()
        self.inpfile.write(" \n")

        materialHeader = """\
$
$ ##########################################################
$                  MATERIAL DEFINITION
$ ##########################################################
$ All material labels present in this model are listed below
$ Need to change the dummy material definition(H2O1) to appropriate one(s) 
$
[MATERIAL]\n"""

        if self.DummyMat:
            self.inpfile.write(materialHeader)
            self.__write_PHITS_source_block__()
            self.inpfile.write(" \n")

        volHeader = """\
$
$ ##########################################################
$                  VOLUME DEFINITION
$ ##########################################################
$ The CAD calculated volume(s) is/are quoted for solid cell(s).
$ Note that the auto-generated void volumes are not calculated,
$ set all at 1.0cm3 tentatively.
$
[VOLUME] off\n"""

        if self.Options["Volume"]:
            self.inpfile.write(volHeader)
            self.__write_PHITS__Volume_block__()
            self.inpfile.write(" \n")

        self.inpfile.close()
        return

    def __write_PHITS_header__(self):

        version = GEOUNED_Version
        releaseDate = GEOUNED_ReleaseDate
        freeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())

        Header = (
            "$ "
            """{}
$   ______ _______  _____      _     _ __   _ _______ ______  
$  |  ____ |______ |     | ___ |     | | \  | |______ |     \ 
$  |_____| |______ |_____|     |_____| |  \_| |______ |_____/
$ Version : {}     {}
$ FreeCAD Version : {} 
$ PHITSFormat Version :  0.0.2.3     06/03/2024\n""".format(
                self.Title, version, releaseDate, freeCAD_Version
            )
        )

        Information = f"""$
$ *************************************************************
$ Original Step file : {self.StepFile}
$
$ Creation Date : {datetime.now()}
$ Solid Cells   : {self.__solidCells__}
$ Total Cells   : {self.__cells__}
$ Surfaces      : {len(self.Surfaces)}
$ Materials     : {len(self.__materials__)}
$
$ **************************************************************


"""

        self.inpfile.write(Header)
        self.inpfile.write(Information)
        return

    def __write_PHITS_cell_block__(self):

        enclenvChk = []
        enclenvChk = self.__stepfile_label_chk__(self.StepFile)

        if enclenvChk:
            print("Unified the inner void cell(s) definition")
            for i, cell in enumerate(self.Cells):
                self.__write_PHITS_cells_uniVoidDef__(cell)
            return
        else:
            for i, cell in enumerate(self.Cells):
                self.__write_PHITS_cells__(cell)
            return

    def __stepfile_label_chk__(self, filename):

        enclenvList = []
        with open(filename) as f:
            enclenvList = f.readlines()
        enclLabel = re.search(
            "enclosure(?P<encl>[0-9]+)_(?P<parent>[0-9]+)_", str(enclenvList)
        )
        envelLabel = re.search(
            "envelope(?P<env>[0-9]+)_(?P<parent>[0-9]+)_", str(enclenvList)
        )
        cond1 = enclLabel == None
        cond2 = envelLabel == None
        return cond1 and cond2

    def __write_PHITS_surface_block__(self):

        for surf in self.Surfaces:
            self.__write_PHITS_surfaces__(surf)

    def __write_PHITS_cells__(self, cell):

        index = cell.label

        # if index is None objet not contain cell definition
        # but a comment to insert between cells
        if cell.__id__ is None:
            comment = self.__commentLine__(cell.Comments)
            self.inpfile.write(comment)
            return

        if cell.Material == 0:
            if cell.MatInfo == "Graveyard":
                cell.MatInfo = "Outer void"
                cellHeader = f"{index:<5d} {-1:<5d}  "

            elif cell.MatInfo == "Graveyard_in":
                cell.MatInfo = "Inner void"
                if self.voidMat != []:
                    self.Materials.add(self.voidMat[0])
                    if abs(self.voidMat[1]) < 1e-2:
                        cellHeader = "{:<5d} {:<5d} {:11.4e} ".format(
                            index, self.voidMat[0], self.voidMat[1]
                        )
                    else:
                        cellHeader = "{:<5d} {:<5d} {:11.7f} ".format(
                            index, self.voidMat[0], self.voidMat[1]
                        )
                else:
                    cellHeader = f"{index:<5d} {0:<5d}  "

            else:
                cellHeader = f"{index:<5d} {0:<5d}  "

        else:
            self.Materials.add(cell.Material)
            if self.Matfile == "" and cell.EnclosureID != 0:
                cellHeader = f"{index:<5d} {cell.Material:<5d} c{cell.Material:<5d} "
            else:
                if abs(cell.Density) < 1e-2:
                    cellHeader = "{:<5d} {:<5d} {:11.4e} ".format(
                        index, cell.Material, cell.Density
                    )
                else:
                    cellHeader = "{:<5d} {:<5d} {:11.7f} ".format(
                        index, cell.Material, cell.Density
                    )

        phitscell = "{}{}\n{}{}".format(
            cellHeader,
            self.__cellFormat__(cell.Definition, offset=len(cellHeader)),
            self.__optionFormat__(cell),
            self.__commentFormat__(cell.Comments, cell.MatInfo),
        )
        self.inpfile.write(phitscell)
        return

    def __write_PHITS_cells_uniVoidDef__(self, cell):

        index = cell.label

        # if index is None objet not contain cell definition
        # but a comment to insert between cells
        if cell.__id__ is None:
            comment = self.__commentLine__(cell.Comments)
            self.inpfile.write(comment)
            return
        """
        # Graveyard_in and Graveyard is changed.
        # Grayveyard sphere is redefined as the boundary surface 
        # of the inner and outer void region.  
        # Auto-generated void cell(s) is/are eliminated.
        # To exclude solid cell(s) from the inner void's defenition,
        # a string of '#(Solid Cell No.)', or inclsolidCells, 
        # is appended to the new inner void cell definition
        # after self.__cellFormat__(cell.Definition) process.       
        """

        if cell.Void:
            if cell.MatInfo == "Graveyard_in":
                cell.MatInfo = "Inner void"

                newInnerVoidCell = cell.Definition.elements[:1]
                cell.Definition.operator = "AND"
                cell.Definition.elements = newInnerVoidCell

                inclSolidCells = ""
                startVoidIndex = self.__solidCells__ + self.startCell
                eliminated_endVoidIndex = self.__cells__ + self.startCell - 3

                if self.startCell == startVoidIndex - 1:
                    inclSolidCells = f"{'':1s}#{self.startCell}"
                else:
                    for i in range(self.startCell, startVoidIndex):
                        inclSolidCells += f"{'':1s}#{i}"

                if startVoidIndex == eliminated_endVoidIndex:
                    one_mervoid_str = "VOID CELL {} merged, so the auto-genarated void definition is eliminated\n"
                    cell.Comments = one_mervoid_str.format(startVoidIndex)
                else:
                    some_mervoid_str = "VOID CELLs {}-{} merged, so the auto-genarated void definitions are eliminated\n"
                    cell.Comments = some_mervoid_str.format(
                        startVoidIndex, eliminated_endVoidIndex
                    )
                if self.voidMat != []:
                    self.Materials.add(self.voidMat[0])
                    if abs(self.voidMat[1]) < 1e-2:
                        cellHeader = "{:<5d} {:<5d} {:11.4e} ".format(
                            index, self.voidMat[0], self.voidMat[1]
                        )
                    else:
                        cellHeader = "{:<5d} {:<5d} {:11.7f} ".format(
                            index, self.voidMat[0], self.voidMat[1]
                        )
                else:
                    cellHeader = f"{index:<5d} {0:<5d}  "

                phitscell = "{}{}\n{}{}\n".format(
                    cellHeader,
                    self.__new_InnerVoid_Def__(
                        inclSolidCells, cell.Definition, offset=len(cellHeader)
                    ),
                    self.__optionFormat__(cell),
                    self.__commentFormat__(cell.Comments, cell.MatInfo),
                )
                self.inpfile.write(phitscell)
                return

            elif cell.MatInfo == "Graveyard":
                cellHeader = f"{index:<5d} {-1:<5d}  "
                cell.MatInfo = "Outer void"

            else:
                return

            """
        # To check auto-generated voids, apply this commented out section instead
        # and comment out above from "if cell.Void:..." to "... else: return"
        # In addition, if you set volCARD = True and want for all void regions to come apperes in [VOLUME],
        # comment out some part in the def __write_PHITS__Volume_block__() section also.       
        if cell.Material == 0:
            print(cell.IsEnclosure)
            if cell.MatInfo == 'Graveyard':
                cellHeader = '{:<5d} {:<5d}  '.format(index,-1)
            else:
                cellHeader = '{:<5d} {:<5d}  '.format(index,0)
           """

        else:
            self.Materials.add(cell.Material)
            if self.Matfile == "" and cell.EnclosureID != 0:
                cellHeader = f"{index:<5d} {cell.Material:<5d} c{cell.Material:<5d} "
            else:
                if abs(cell.Density) < 1e-2:
                    cellHeader = "{:<5d} {:<5d} {:11.4e} ".format(
                        index, cell.Material, cell.Density
                    )
                else:
                    cellHeader = "{:<5d} {:<5d} {:11.7f} ".format(
                        index, cell.Material, cell.Density
                    )

        phitscell = "{}{}\n{}{}".format(
            cellHeader,
            self.__cellFormat__(cell.Definition, offset=len(cellHeader)),
            self.__optionFormat__(cell),
            self.__commentFormat__(cell.Comments, cell.MatInfo),
        )
        self.inpfile.write(phitscell)
        return

    def __write_PHITS_surfaces__(self, surface):
        """Write the surfaces in PHITS format"""

        PHITS_def = PHITSSurface(surface.Index, surface.Type, surface.Surf)
        if PHITS_def:
            PHITS_def += "\n"
            self.inpfile.write(PHITS_def)
        else:
            print(f"Surface {surface.Type} cannot be written in PHITS input")
        return

    def __write_PHITS_source_block__(self):

        if self.DummyMat:
            mat = list(self.Materials)
            mat.sort()
            MATID = []
            MATCARD = ""
            for i, cell in enumerate(self.Cells):
                if (cell.Material in mat) and (cell.Material not in MATID):
                    MATID.append(cell.Material)
                    if self.Matfile == "" and cell.EnclosureID != 0:
                        mismat_comment = "$ Change dummyMat M{}, {} c{} g/cm3 is assigned\n M{:<6d} H 2 O 1\n"
                        MATCARD += mismat_comment.format(
                            cell.Material, cell.MatInfo, cell.Material, cell.Material
                        )
                    else:
                        mat_comment = "$ Change dummyMat M{} to {}, Density = {}g/cm3\n M{:<6d} H 2 O 1\n"
                        MATCARD += mat_comment.format(
                            cell.Material,
                            cell.MatInfo,
                            -1 * cell.Density,
                            cell.Material,
                        )
            Block = MATCARD + "\n"

        self.inpfile.write(Block)

    def __write_PHITS__Volume_block__(self):

        vol = f"{'':5s}reg{'':5s}vol\n"

        startVoidIndex = self.__solidCells__ + self.startCell
        eliminated_endVoidIndex = self.__cells__ + self.startCell - 3

        enclenvChk = []
        enclenvChk = self.__stepfile_label_chk__(self.StepFile)

        if enclenvChk and self.Options["Volume"]:
            for i, cell in enumerate(self.Cells):
                if cell.__id__ is not None:
                    if cell.Void and startVoidIndex == eliminated_endVoidIndex:
                        if cell.label == startVoidIndex:
                            print(
                                "Eliminated the merged void cell {} from [VOLUME] section".format(
                                    cell.label
                                )
                            )
                        else:
                            vol += f"{'':6s}{cell.label}{'':6s}1.0\n"
                    elif cell.Void:
                        if cell.label in range(
                            startVoidIndex, eliminated_endVoidIndex + 1
                        ):
                            print(
                                "Eliminated the merged void cell {} from [VOLUME] section".format(
                                    cell.label
                                )
                            )
                        else:
                            vol += f"{'':6s}{cell.label}{'':6s}1.0\n"
                    else:
                        vol += f"{'':6s}{cell.label}{'':6s}{cell.Volume * 0.001:6e}\n"
        else:
            if self.Options["Volume"]:
                for i, cell in enumerate(self.Cells):
                    if cell.__id__ is not None:
                        if cell.Void:
                            vol += f"{'':6s}{cell.label}{'':6s}1.0\n"
                        else:
                            vol += "{:6s}{}{:6s}{:6e}\n".format(
                                "", cell.label, "", cell.Volume * 1e-3
                            )

        self.inpfile.write(vol)

        """
        # If you want for all eliminated the merged void cells to come apperes in [VOLUME],
        # apply this commented out section instead
        # and comment out above from "startVoidIndex =..." to "... self.inpfile.write(vol)"
        if self.Options['Volume']:
            for i,cell in enumerate(self.Cells):
                if cell.__id__ is not None :
                    if cell.Void : 
                        vol +='{:6s}{}{:6s}1.0\n'.format('',cell.label,'')
                    else:
                        vol +='{:6s}{}{:6s}{:e}\n'.format('',cell.label,'',cell.Volume*1e-3)

        self.inpfile.write(vol)
        """

    def __cellFormat__(self, Definition, offset=11):
        return writePHITSCellDef(Definition, tabspace=11, offset=offset)

    def __new_InnerVoid_Def__(self, innerSolidCells, Definition, offset=11):
        newInnerVoidDef = self.__cellFormat__(Definition, offset)
        strdef = CellString(tabspace=11)
        strdef.add(newInnerVoidDef + innerSolidCells)
        strdef.wrapLine(offset)
        return strdef.str

    def __optionFormat__(self, cell):

        option = ""
        if self.Options["Volume"]:
            if not cell.Void:
                option = f"${'':11s}Vol={cell.Volume * 0.001:e} cm3\n"
            else:
                option = f"${'':11s}Vol=1.0 cm3\n"

        if self.Options["Universe"] is not None:
            option += f"{'':11s}U={self.Options['Universe']}\n"

        return option

    def __commentFormat__(self, cComment, mComment=None):

        comment = ""
        if mComment:
            mComment = mComment.split("\n")
            for c in mComment:
                if c:
                    comment += f"{'':11s}${c}\n"

        if cComment.strip() != "":
            cComment = cComment.strip().split("\n")
            for c in cComment:
                if c:
                    comment += f"{'':11s}${c}\n"
        return comment

    def __commentLine__(self, lineComment):
        lineComment = lineComment.strip().split("\n")
        comment = ""
        if lineComment:
            comment = "$ \n"
            for c in lineComment:
                if c:
                    comment += f"$ {c}\n"
            comment += "$ \n"
        return comment

    def __getSurfaceTable__(self):
        self.surfaceTable = {}
        self.__solidCells__ = 0
        self.__cells__ = 0
        self.__materials__ = set()

        for i, CellObj in enumerate(self.Cells):
            if CellObj.__id__ is None:
                continue
            self.__cells__ += 1
            if CellObj.Material != 0:
                self.__materials__.add(CellObj.Material)

            surf = CellObj.Definition.getSurfacesNumbers()
            if not CellObj.Void:
                self.__solidCells__ += 1
            for index in surf:
                if index in self.surfaceTable.keys():
                    self.surfaceTable[index].add(i)
                else:
                    self.surfaceTable[index] = {i}
        return

    def __simplifyPlanes__(self, Surfaces):

        for p in Surfaces["PX"]:
            if p.Surf.Axis[0] < 0:
                p.Surf.Axis = FreeCAD.Vector(1, 0, 0)
                self.__changeSurfSign__(p)

        for p in Surfaces["PY"]:
            if p.Surf.Axis[1] < 0:
                p.Surf.Axis = FreeCAD.Vector(0, 1, 0)
                self.__changeSurfSign__(p)

        for p in Surfaces["PZ"]:
            if p.Surf.Axis[2] < 0:
                p.Surf.Axis = FreeCAD.Vector(0, 0, 1)
                self.__changeSurfSign__(p)

        if opt.prnt3PPlane:
            for p in Surfaces["P"]:
                if p.Surf.pointDef:
                    axis, d = pointsToCoeffs(p.Surf.Points)
                    if isOposite(axis, p.Surf.Axis):
                        self.__changeSurfSign__(p)
        return

    def __sortedSurfaces__(self, Surfaces):
        temp = SurfacesDict(Surfaces)
        surfList = []
        for ind in range(
            Surfaces.IndexOffset, Surfaces.surfaceNumber + Surfaces.IndexOffset
        ):
            s = temp.getSurface(ind + 1)
            if s is not None:
                surfList.append(s)
                temp.delSurface(ind + 1)
        return surfList

    def __changeSurfSign__(self, p):

        if p.Index not in self.surfaceTable.keys():
            print(
                f"{p.Type} Surface {p.Index} not used in cell definition)",
                p.Surf.Axis,
                p.Surf.Position,
            )
            return

        for ic in self.surfaceTable[p.Index]:
            surf = self.Cells[ic].Definition.getSurfacesNumbers()
            for s in surf:
                if s == p.Index:
                    changeSurfSign(s, self.Cells[ic].Definition)

    def __get_solidCellVolume__(self):

        solidList = []
        volumeList = []
        for m in self.Cells:
            if m.CellType == "solid" and m.__id__ is not None:
                solidList.append(m.label)
                volumeList.append(m.Volume * 1e-3)
        return solidList, volumeList
