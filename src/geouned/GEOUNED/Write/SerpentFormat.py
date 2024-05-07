#################################
# Module to write Serpent input #
#################################
from datetime import datetime

import FreeCAD

from ..CodeVersion import *
from ..Utils.BasicFunctions_part1 import is_opposite, points_to_coeffs
from ..Utils.Functions import SurfacesDict
from .Functions import serpent_surface, change_surf_sign, write_serpent_cell_def


class SerpentInput:
    def __init__(
        self,
        Meta,
        Surfaces,
        title,
        volSDEF,
        volCARD,
        UCARD,
        dummyMat,
        prnt3PPlane,
        StepFile,
        tolerances,
        numeric_format,
        options,
    ):
        self.prnt3PPlane = prnt3PPlane
        self.StepFile = StepFile
        self.tolerances = tolerances
        self.numeric_format = numeric_format
        self.Title = title
        self.VolSDEF = volSDEF
        self.VolCARD = volCARD
        self.U0CARD = UCARD
        self.dummyMat = dummyMat
        self.Cells = Meta
        self.options = options

        # TODO allow user to edit these from CadToCsg class
        self.particle = ("n", "p")
        self.part = "p"

        self.__get_surface_table__()
        self.__simplify_planes__(Surfaces, prnt3PPlane)

        self.Surfaces = self.__sorted_surfaces__(Surfaces)
        self.Materials = set()

        return

    # def setSRC(self,data):

    #     if data[0] is not None:
    #       sdef     = f'src point {self.part}\n'
    #       self.src_sphere = (sdef)
    #     else:
    #       self.src_sphere = None
    #     xmin,xmax,ymin,ymax,zmin,zmax = data[1]

    #     sdef = 'SDEF PAR={} X=D1 Y=D2 Z=D3 \n'.format(self.part)
    #     SI1  = 'SI1 {:13.7e} {:13.7e} \n'.format(xmin*0.1,xmax*0.1)
    #     SI2  = 'SI2 {:13.7e} {:13.7e} \n'.format(ymin*0.1,ymax*0.1)
    #     SI3  = 'SI3 {:13.7e} {:13.7e} \n'.format(zmin*0.1,zmax*0.1)
    #     SP1  = 'SP1 0  1 \n'
    #     SP2  = 'SP2 0  1 \n'
    #     SP3  = 'SP3 0  1 \n'
    #     self.SDEF_box    = (sdef,SI1,SI2,SI3,SP1,SP2,SP3)

    def write_input(self, filename, verbose):
        if verbose:
            print(f"write Serpent file {filename}")
        self.inpfile = open(filename, "w", encoding="utf-8")
        self.__write_header__()
        cellblockHeader = """\
% --- CELL DEFINITIONS 
"""
        self.inpfile.write(cellblockHeader)
        self.__write_cell_block__()
        self.inpfile.write(" \n")

        surfaceHeader = """\
% --- SURFACE DEFINITIONS 
"""
        self.inpfile.write(surfaceHeader)
        self.__write_surface_block__()
        self.inpfile.write(" \n")

        self.__write_source_block__()

        self.inpfile.close()
        return

    def __write_header__(self):

        version = GEOUNED_Version
        releaseDate = GEOUNED_ReleaseDate
        freeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())

        Header = f"""{self.Title}
%   ______ _______  _____      _     _ __   _ _______ ______  
%  |  ____ |______ |     | ___ |     | | \\  | |______ |     \\ 
%  |_____| |______ |_____|     |_____| |  \\_| |______ |_____/
% Version : {version}     {releaseDate}
% FreeCAD Version : {freeCAD_Version} 
"""

        Information = f"""%
% *************************************************************
% Original Step file : {self.StepFile}
%
% Creation Date : {datetime.now()}
% Solid Cells   : {self.__solidCells__}
% Total Cells   : {self.__cells__}
% Surfaces      : {len(self.Surfaces)}
% Materials     : {len(self.__materials__)}
%
% **************************************************************
"""
        self.inpfile.write(Header)
        self.inpfile.write(Information)
        return

    def __write_surface_block__(self):

        for surf in self.Surfaces:
            self.__write_surfaces__(surf)

    def __write_cell_block__(self):

        for i, cell in enumerate(self.Cells):
            self.__write_cells__(cell)
        return

    # to_mod
    def __write_cells__(self, cell):
        index = cell.label

        # If index is None, the object does not contain cell definition
        # but a comment to insert between cells
        if cell.__id__ is None:
            comment = self.comment_line(cell.Comments)
            self.inpfile.write(comment)
            return

        if self.U0CARD is not None:
            if cell.Material == 0:
                cellHeader = f'cell {index:<5d} {self.U0CARD} {"void":<5d}  '
            else:
                self.Materials.add(cell.Material)
                cellHeader = f"cell {index:<5d} {self.U0CARD} {cell.Material:<5d} "

            serpent_cell = (
                f"{cellHeader}{self.__cell_format__(cell.Definition, offset=len(cellHeader))}"
                f"{self.comment_format(cell.Comments, cell.MatInfo)}"
            )
            self.inpfile.write(serpent_cell)
        else:
            if cell.Material == 0 and cell.MatInfo != "Graveyard":
                cellHeader = f"cell {index:<5d} 0 void "
            elif cell.Material == 0 and cell.MatInfo == "Graveyard":
                cellHeader = f"cell {index:<5d} 0 outside "
            else:
                self.Materials.add(cell.Material)
                cellHeader = f"cell {index:<5d}  0  {cell.Material:<5d} "

        serpent_cell = (
            f"{cellHeader}{self.__cell_format__(cell.Definition, offset=len(cellHeader))}"
            f"{self.comment_format(cell.Comments, cell.MatInfo)}"
        )
        self.inpfile.write(serpent_cell)

        return

    def __write_surfaces__(self, surface):
        """Write the surfaces in Serpent format"""

        Serpent_def = serpent_surface(
            surface.Index,
            surface.Type,
            surface.Surf,
            self.tolerances,
            self.numeric_format,
            self.options,
        )
        if Serpent_def:
            Serpent_def += "\n"
            self.inpfile.write(Serpent_def)
        else:
            print(f"Surface {surface.Type} cannot be written in Serpent input")
        return

    # No void all option in Serpent. For now remove addition of source.

    def __write_source_block__(self):

        #       if self.SDEF_sphere is None:  return
        MODE = f"\nset nps 1e6\nset bc 1"
        if self.dummyMat:
            mat = list(self.Materials)
            mat.sort()
            MATCARD = ""
            for m in mat:
                MATCARD += f"mat {m:<6d} {self.cell.Density:11.4e} \n1001 1 \n"
            Block = MATCARD + "% \n" + MODE
        else:
            Block = MODE

        # Not included for Serpent
        #    if self.VolSDEF:
        #       Block += 'PRDMP 2J -1\n'
        #       for line in self.SDEF_box :
        #          Block += 'C ' + line
        #       for line in self.SDEF_sphere:
        #          Block += line

        #       celList,volList = self.__get_solid_cell_volume__()

        #       F4Tally = CardLine('F4:{} '.format(self.part))
        #       F4Tally.extend(celList)
        #       SD4     = CardLine('SD4  ',fmt='13.7e')
        #       SD4.extend(volList)

        #       Block += F4Tally.str+'\n'
        #       if not self.VolCARD :
        #          Block += SD4.str
        #       else :
        #          Block += 'C Cell volume normalization is set in cell cards VOL\n'
        #    else :
        #       for line in self.SDEF_sphere :
        #          Block += 'C ' + line
        #       for line in self.SDEF_box :
        #          Block +=  line

        self.inpfile.write(Block)

    def __cell_format__(self, Definition, offset=11):
        return write_serpent_cell_def(Definition, tabspace=11, offset=offset)

    # Function not relevant for Serpent : No importance setting, universes assigned elsewhere.
    # Volumes only defined on tally cards.
    # def __option_format__(self,cell):

    #    option = ''
    #    if self.VolCARD:
    #       if not cell.Void :
    #         option ='{:11s}Vol={:e}\n'.format('',cell.Volume*1e-3)
    #       else:
    #         option = '{:11s}Vol=1.0\n'.format('')

    #    option +=  '{:11s}'.format('')
    #    for p in self.Options['Particle']:
    #      if cell.MatInfo == 'Graveyard' :
    #        option +=  'imp:{}=0     '.format(p)
    #      else:
    #        option +=  'imp:{}=1.0   '.format(p)

    #    if self.Options['Universe'] is not None:
    #         option += 'U={}  '.format(self.Options['Universe'])
    #    option += '\n'

    #    return option

    def comment_format(self, cComment, mComment=None):
        comment = ""
        if mComment:
            mComment = mComment.split("\n")
            for c in mComment:
                if c:
                    comment += f"{'':11s}%{c}\n"

        if cComment.strip() != "":
            cComment = cComment.strip().split("\n")
            for c in cComment:
                if c:
                    comment += f"{'':11s}%{c}\n"
        return comment

    def comment_line(self, lineComment):
        lineComment = lineComment.strip().split("\n")
        comment = ""
        if lineComment:
            comment = "% \n"
            for c in lineComment:
                if c:
                    comment += f"% {c}\n"
            comment += "% \n"
        return comment

    def __get_surface_table__(self):
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

            surf = CellObj.Definition.get_surfaces_numbers()
            if not CellObj.Void:
                self.__solidCells__ += 1
            for index in surf:
                if index in self.surfaceTable.keys():
                    self.surfaceTable[index].add(i)
                else:
                    self.surfaceTable[index] = {i}
        return

    def __simplify_planes__(self, Surfaces, prnt3PPlane):

        for p in Surfaces["PX"]:
            if p.Surf.Axis[0] < 0:
                p.Surf.Axis = FreeCAD.Vector(1, 0, 0)
                self.__change_surf_sign__(p)

        for p in Surfaces["PY"]:
            if p.Surf.Axis[1] < 0:
                p.Surf.Axis = FreeCAD.Vector(0, 1, 0)
                self.__change_surf_sign__(p)

        for p in Surfaces["PZ"]:
            if p.Surf.Axis[2] < 0:
                p.Surf.Axis = FreeCAD.Vector(0, 0, 1)
                self.__change_surf_sign__(p)

        if prnt3PPlane:
            for p in Surfaces["P"]:
                if p.Surf.pointDef:
                    axis, d = points_to_coeffs(p.Surf.Points)
                    if is_opposite(axis, p.Surf.Axis):
                        self.__change_surf_sign__(p)

        return

    def __sorted_surfaces__(self, Surfaces):
        temp = SurfacesDict(Surfaces)
        surfList = []
        for ind in range(
            Surfaces.IndexOffset, Surfaces.surfaceNumber + Surfaces.IndexOffset
        ):
            s = temp.get_surface(ind + 1)
            if s is not None:
                surfList.append(s)
                temp.del_surface(ind + 1)
        return surfList

    def __change_surf_sign__(self, p):

        if p.Index not in self.surfaceTable.keys():
            print(
                f"{p.Type} Surface {p.Index} not used in cell definition)",
                p.Surf.Axis,
                p.Surf.Position,
            )
            return

        for ic in self.surfaceTable[p.Index]:
            surf = self.Cells[ic].Definition.get_surfaces_numbers()
            for s in surf:
                if s == p.Index:
                    change_surf_sign(s, self.Cells[ic].Definition)

    def __get_solid_cell_volume__(self):

        solidList = []
        volumeList = []
        for m in self.Cells:
            if m.CellType == "solid" and m.__id__ is not None:
                solidList.append(m.label)
                volumeList.append(m.Volume * 1e-3)

        return solidList, volumeList
