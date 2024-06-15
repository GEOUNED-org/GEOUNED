import typing
import FreeCAD
import Import

from pathlib import Path

from .Modules.buildCAD import buildCAD, makeTree
from .Modules.MCNPinput import McnpInput
from .Modules.Objects import CadCell
from .Modules.XMLinput import XmlInput


class CsgToCad:
    """Base class for the conversion of CSG to CAD models"""

    def __init__(self):
        pass

    def export_cad(
        self,
        input_filename: str,
        csg_format: str,
        output_filename: str = "cad_from_csg",
        bounding_box: typing.Tuple[int, int, int, int, int, int] = (-1000, -1000, -1000, 1000, 1000, 1000),
        universe_start: int = 0,
        level_max: str = "all",
        cell_range_type: str = "all",
        cell_range: typing.Tuple[int] = (),
        mat_range_type: str = "all",
        mat_range: typing.Tuple[int] = (),
        # TODO add these to the method signature
        # splitTolerance in the Options class
        # mat =  this is in the CADselection dictionary but not in the docs https://github.com/GEOUNED-org/GEOUNED/blob/76ef697c7dca6a828c7498996ff3313859c872f2/docs/User_Guide_GEOUNED_v2.0.pdf
    ):
        """export the CSG geometry in OpenMC or MCNP format to a CAD model.

        Args:
            input_filename (str): The filename and path of the input CSG text file.
            csg_format (str): The format of the CSG input file, options are 'openmc_xml' or 'mcnp'
            output_filename (str, optional): The filename stem and path of the output file created.
                Two files will be created with the '.step' suffix and one with the 'FCStd' suffix.
                Defaults to 'cad_from_csg'.
            bounding_box (typing.Tuple[int, int, int, int, int, int], optional): The bounding box
                coordinates of the CSG geometry. This should encompass the entire CSG geometry.
                Format is (xmin, ymin, zmin, xmax, ymax, zmax) Defaults to (-1000, -1000, -1000,
                1000, 1000, 1000).
            universe_start (int, optional): The Universe ID to be converted to CAD. If universe_start
                is left as 0 then all the universes and any nested universes are converted. If set
                then the universe and all its nested universes are converted. Defaults to 0.
            level_max (str, optional): Level maximum of nested universe to be translated. If
                level_max < highest nested universe level, cells inside the container cell whose
                level is level_max will not be translated. This container cell will be the CAD
                solid written in the CAD file. Defaults to "all".
            cell_range_type (str, optional): Define how to consider the range values.
                setting to 'all' results in all the cells with any cell ID being converted (range a no effect).
                setting to 'include' results in only cells defined in 'range' being converted.
                setting to 'exclude' results in exclude all cells defined in range. Defaults to "all".
            cell_range (typing.Tuple[int], optional): list of cells to be included/excluded for the conversion.
                Defaults to ().
            mat_range_type (str, optional): Define how to consider the range values.
                setting to 'all' results in all the materials with any cell ID being converted (range a no effect).
                setting to 'include' results in only materials defined in 'range' being converted.
                setting to 'exclude' results in exclude all materials defined in range. Defaults to "all".
            mat_range (typing.Tuple[int], optional): list of materials to be included/excluded for the conversion.
                Defaults to ().

        Raises:
            ValueError: If the csg_format is not 'openmc_xml' or 'mcnp' then a ValueError is raised.
        """

        # TODO check file extensions are correct
        # if Path(output_filename).suffix not in ['.stp', '.step']:
        #     raise ValueError(f"output file must have a .stp or .step extension, not {universe_start.suffix}")

        # get geometry definition from OpenMC XML or MCNP input
        if csg_format == "mcnp":
            geo = McnpInput(input_filename)
        elif csg_format == "openmc_xml":
            geo = XmlInput(input_filename)
        else:
            msg = f"input format type {csg_format} is not supported. Supported options are 'openmc_xml' or 'mcnp'"
            raise ValueError(msg)

        Path(output_filename).parent.mkdir(parents=True, exist_ok=True)

        UnivCell = CadCell()
        UnivCell.shape = UnivCell.makeBox(FreeCAD.BoundBox(*bounding_box))

        # TODO make these variable names lower case in the downstream code

        CADselection = {
            "Ustart": universe_start,
            "levelMax": level_max,
            "cell": [cell_range_type, cell_range],
            "mat": [mat_range_type, mat_range],
            "format": csg_format,
            "cell_range_type": cell_range_type,
            "cell_range": cell_range,
            "mat_range_type": mat_range_type,
            "mat_range": mat_range,
        }

        # TODO don't return fails variable, just fail in the method and raise the error there
        CADCells, fails = buildCAD(UnivCell, geo, CADselection)

        if fails:
            print("failed in conversion", fails)

        CADdoc = FreeCAD.newDocument("converted_with_geouned")

        makeTree(CADdoc, CADCells)
        Import.export(CADdoc.Objects[0:1], f"{output_filename}.step")
        CADdoc.saveAs(f"{output_filename}.FCStd")
