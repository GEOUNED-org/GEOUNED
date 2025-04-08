import FreeCAD
import Import
import typing

from pathlib import Path

from .Modules.buildCAD import makeTree, AssignSurfaceToCell, BuildUniverseCells
from .Modules.Utils.booleanFunction import BoolSequence
from .Modules.Utils.boundBox import solid_plane_box
from .Modules.data_class import BoxSettings
from .Modules.Objects import CadCell
from .Modules.MCNPinput import McnpInput
from .Modules.XMLinput import XmlInput


class CsgToCad:
    """Base class for the conversion of CSG to CAD models

    Args:
        BoxSettings (geouned.BoxSettings, optional): Adjust the default parameters
        for the solid boundbox creation. Defaults to a geouned.BoxSettings with default
        attributes values.
    """

    def __init__(self, settings: BoxSettings = BoxSettings()):
        self.settings = settings
        self.cell_range_type = "all"
        self.cell_range = None
        self.mat_range_type = "all"
        self.mat_range = None
        self.buildCAD_list = []

    def read_csg_file(self, input_filename: str, csg_format: str):
        """Reads the geometry definition from MCNP or OpenMC XML input.

        Args:
            input_filename (str): The filename and path of the input CSG text file.
            csg_format (str): The format of the CSG input file, options are 'mcnp' or 'openmc_xml'

        Raises:
            ValueError: If the csg_format is not 'openmc_xml' or 'mcnp' then a ValueError is raised."""

        if csg_format == "mcnp":
            self.geometry = McnpInput(input_filename)
        elif csg_format == "openmc_xml":
            self.geometry = XmlInput(input_filename)
        else:
            msg = f"input format type {csg_format} is not supported. Supported options are 'openmc_xml' or 'mcnp'"
            raise ValueError(msg)
            # read all surfaces definition

        self.input_filename = input_filename
        self.geometry.GetSurfaces()  # scale units change are carried out in GetSurfaces method
        self.geometry.GetLevelStructure()

    def cell_filter(self, type: str = "all", cells: typing.Union[None, list, tuple] = None):
        """Selects the cells to build from the CSG geometry in MCNP or OpenMC format and export to the CAD model.

        Args:
            type (str, optional): Filtering type. Allowed values "all", "include", "exclude". Default to all.
            cells (None, list, tuple, optional): List of cells to include or exclude. If type is "all" has no effect. Default to None
        """

        if type in ("exclude", "include", "all"):
            self.cell_range_type = type
        else:
            self.cell_range_type = "all"
            print("bad cells range type. Ignored")

        if cells is not None:
            self.cell_range = cells[:]
        else:
            if self.cell_range_type == "exclude":
                self.cell_range_type = "all"

    def material_filter(self, type="all", materials=None):
        """Selects the materials of the cells to build from the CSG geometry in MCNP or OpenMC format and export to the CAD model.

        Args:
            type (str, optional): Filtering type. Allowed values "all", "include", "exclude". Default to all.
            materials (None, list, tuple, optional): List of cells to include or exclude. If type is "all" has no effect. Default to None
        """

        if type in ("exclude", "include", "all"):
            self.mat_range_type = type
        else:
            self.mat_range_type = "all"
            print("bad cells range type. Ignored")

        if materials is not None:
            self.mat_range = materials[:]
        else:
            if self.mat_range_type == "exclude":
                self.mat_range_type = "all"

    def build_container(self, cell_label: int, depth: int = -1):
        """Build the universe contained in the cell "cell_label". The level of nested universes to consider is
        controlled by the parameter "depth". The universe is located inside the container cell according to the transformation.

        Args:
            cell_label (int): Label of the cell containing the universe to convert.
            depth (int, optional): Depth level of the nested to considered. Default to -1 (all nested universes)."""

        print(f"Build container cell {cell_label}")
        # get and build container universe
        UnivCell = self.geometry.GetCell(cell_label, self.settings)
        UnivCell.definition = BoolSequence(UnivCell.definition.str)

        UnivCell.build_BoundBox(enlarge=0.2)
        if UnivCell.boundBox.Orientation == "Forward" and UnivCell.boundBox.Box is None:
            UnivCell.shape = None
            print(f"Cell {UnivCell.name} BoundBox is null")
            return
        else:
            if UnivCell.boundBox.Orientation == "Forward":
                UnivCell.externalBox = UnivCell.boundBox

        debug = True
        if debug:
            UnivCell.buildShape(simplify=False)
        else:
            try:
                UnivCell.buildShape(simplify=False)
            except:
                print(f"fail converting cell {UnivCell.name}")

        # generate universe cells inside container
        matcel_list = {
            "mat": (self.mat_range_type, self.mat_range),
            "cell": (self.cell_range_type, self.cell_range),
        }

        UniverseCells, modelSurfaces = self.geometry.GetFilteredCells(UnivCell.FILL, depth, matcel_list, self.settings)
        AssignSurfaceToCell(UniverseCells, modelSurfaces)

        Ustart = UnivCell.FILL
        UnivCell.level = None
        for lev, Univ in self.geometry.levels.items():
            if Ustart in Univ:
                UnivCell.level = lev
                break

        if depth == -1:
            levelMax = len(self.geometry.levels)
        else:
            levelMax = min(lev + depth, len(self.geometry.levels))

        startInfo = (Ustart, levelMax)
        CADCells, fails = BuildUniverseCells(startInfo, UnivCell, UniverseCells, universeCut=True)
        if fails:
            print("failed cell conversion:", fails)
        self.buildCAD_list.append(CADCells)

    def build_universe(self, U: typing.Union[None, int] = None, depth: int = -1):
        """Build the universe U. The level of nested universes is to consider is controlled by the parameter "depth".
        The universe is build in its own coordinate system.

        Args:
            U (int, optional): Value of the universe to convert. Default to None (root universe).
            depth (int, optional): Depth level of the nested to considered. Default to -1 (all nested universes)."""

        UniverseCut = True
        # UnivCell is the virtual(infinite) cell container
        UnivCell = CadCell(settings=self.settings)
        if U == None:
            root_universe = self.geometry.levels[0][0]
        else:
            root_universe = U

        UnivCell.FILL = root_universe
        UnivCell.name = None
        UnivCell.MAT = None

        # read Cells and group into universes
        matcel_list = {
            "mat": (self.mat_range_type, self.mat_range),
            "cell": (self.cell_range_type, self.cell_range),
        }
        # get all cells of the universe U and subuniverses in the FILL cell of universe U
        # depth level for searching for nested universes
        UniverseCells, modelSurfaces = self.geometry.GetFilteredCells(root_universe, depth, matcel_list, self.settings)

        # assign to each cell the surfaces belonging to the cell
        AssignSurfaceToCell(UniverseCells, modelSurfaces)

        Ustart = root_universe
        UnivCell.level = None
        for lev, Univ in self.geometry.levels.items():
            if Ustart in Univ:
                UnivCell.level = lev
                break

        if depth == -1:
            levelMax = len(self.geometry.levels) - 1
        else:
            levelMax = min(lev + depth, len(self.geometry.levels) - 1)

        startInfo = (Ustart, levelMax)
        CADCells, fails = BuildUniverseCells(startInfo, UnivCell, UniverseCells, universeCut=UniverseCut)
        self.buildCAD_list.append(CADCells)
        if fails:
            print("failed cell conversion:", fails)

    def export_cad(self, output_filename: str = ""):
        """export the CSG geometry in OpenMC or MCNP format to a CAD model.

        Args:
            output_filename (str, optional): The filename stem and path of the output file created.
                Two files will be created with the '.step' suffix and one with the 'FCStd' suffix.
                Defaults to name of the csg file + stp.
        """

        if output_filename == "":
            output_filename = Path(self.input_filename).name

        Path(output_filename).parent.mkdir(parents=True, exist_ok=True)

        fullname = Path(output_filename).name
        suffix = Path(output_filename).suffix
        if suffix != "":
            barename = fullname[0 : -len(suffix)]
            output_filename = output_filename[0 : -len(suffix)]
        else:
            barename = fullname

        if suffix not in (".stp", ".step"):
            suffix = ".stp"

        CADdoc = FreeCAD.newDocument("converted_with_geouned")

        CADobj = CADdoc.addObject("App::Part", "Universes")
        CADobj.Label = barename

        for CAD in self.buildCAD_list:
            CADobj.addObject(makeTree(CADdoc, CAD))

        Import.export(CADdoc.Objects[0:1], output_filename + suffix)
        CADdoc.saveAs(f"{output_filename}.FCStd")
