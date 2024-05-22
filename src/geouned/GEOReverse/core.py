import typing
import FreeCAD
import Import

from pathlib import Path

from .Modules.buildCAD import buildCAD, makeTree
from .Modules.MCNPinput import McnpInput
from .Modules.Objects import CadCell
from .Modules.XMLinput import XmlInput


def CsgToCad():
    
    def __init__(
        self,
    ):
        pass

    def export_cad(
        self,
        input_filename,
        output_filename,
        csg_format: str,
        out_box: typing.Tuple[int, int, int, int, int, int] = (-100, -100, -100, 100, 100, 100),
        universe_start: int = 0,
        level_max: str = "all",
        # rangeType
        # range
        # splitTolerance in the Options class
        # mat =  this is in the CADselection dictionary but not in the docs https://github.com/GEOUNED-org/GEOUNED/blob/76ef697c7dca6a828c7498996ff3313859c872f2/docs/User_Guide_GEOUNED_v2.0.pdf
    ):

        if Path(output_filename).suffix not in ['.stp', '.step']:
            raise ValueError(f"output file must have a .stp or .step extension, not {universe_start.suffix}")

        # get geometry definition from OpenMC XML or MCNP input
        if csg_format == "mcnp":
            geo = McnpInput(input_filename)
        elif csg_format == "openmc_xml":
            geo = XmlInput(input_filename)
        else:
            msg = f"input format type {csg_format} is not supported. Supported options are 'openmc_xml' or 'mcnp'"
            raise ValueError(msg)
    
        UnivCell = CadCell()
        UnivCell.shape = UnivCell.makeBox(FreeCAD.BoundBox(*out_box))

        # TODO make these variable names lower case in the downstream code
        CADselection = {
            "Ustart": universe_start,
            "levelMax": level_max,
            # "cell": cell
            # # "mat": mat
            "format": csg_format,
        }

        # TODO don't return fails varible, just fail in the method and raise the error there
        CADCells, fails = buildCAD(UnivCell, geo, CADselection)

        if fails:
            print("failed in conversion", fails)

        CADdoc = FreeCAD.newDocument("converted_with_geouned")

        makeTree(CADdoc, CADCells)
        Import.export(CADdoc.Objects[0:1], universe_start)
        CADdoc.saveAs(universe_start + ".FCStd")