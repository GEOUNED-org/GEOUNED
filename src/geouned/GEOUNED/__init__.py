# Init file of GEOUNED module
#

# We load the STEP and the materials

import configparser
import logging
import typing
import json
from datetime import datetime
from os import mkdir, path
from typing import get_type_hints

import FreeCAD
import Part
from tqdm import tqdm
from pathlib import Path

from .CodeVersion import *
from .Conversion import CellDefinition as Conv
from .Cuboid.translate import translate
from .Decompose import Decom_one as Decom
from .LoadFile import LoadSTEP as Load
from .Utils import Functions as UF
from .Utils.BooleanSolids import build_c_table_from_solids
from .Void import Void as Void
from .Write.Functions import write_mcnp_cell_def
from .Write.WriteFiles import write_geometry
from .Utils.data_classes import Options, Tolerances, NumericFormat, Settings

logger = logging.getLogger(__name__)


class CadToCsg:
    """Base class for the conversion of CAD to CSG models

    Args:
        stepFile (str, optional): Name of the CAD file (in STEP format) to
            be converted. Defaults to "".
        options (geouned.Options, optional): An instance of a geouned.Options
            class with the attributes set for the desired conversion. Defaults
            to a geouned.Options with default attributes values.
        tolerances (geouned.Tolerances, optional): An instance of a
            geouned.Tolerances class with the attributes set for the desired
            conversion. Defaults to a geouned.Tolerances with default
            attributes values.
        numeric_format (geouned.NumericFormat, optional): An instance of a
            geouned.NumericFormat class with the attributes set for the desired
            conversion. Defaults to a geouned.NumericFormat with default
            attributes values.
        settings (geouned.Settings, optional): An instance of a
            geouned.Settings class with the attributes set for the desired
            conversion. Defaults to a geouned.Settings with default
            attributes values.
    """

    def __init__(
        self,
        stepFile: str = "",
        options: Options = Options(),
        tolerances: Tolerances = Tolerances(),
        numeric_format: NumericFormat = NumericFormat(),
        settings: Settings = Settings(),
    ):

        self.stepFile = stepFile
        self.options = options
        self.tolerances = tolerances
        self.numeric_format = numeric_format
        self.settings = settings

    def export_csg(
        self,
        title: str = "Converted with GEOUNED",
        geometryName: str = "csg",
        outFormat: typing.Tuple[str] = (
            "openMC_XML",
            "openMC_PY",
            "serpent",
            "phits",
            "mcnp",
        ),
        volSDEF: bool = False,
        volCARD: bool = True,
        UCARD: int = 101,
        dummyMat: bool = False,
        cellCommentFile: bool = False,
        cellSummaryFile: bool = True,
    ):
        """Writes out a CSG file in the requested Monte Carlo code format.

        Args:
            title (str, optional): Title of the model written at the top of the
                output file. Defaults to "Geouned conversion".
            geometryName (str, optional): the file stem of the output file(s).
                Defaults to "converted_with_geouned".
            outFormat (typing.Tuple[str], optional): Format for the output
                geometry. Available format are: "mcnp", "openMC_XML",
                "openMC_PY", "phits" and "serpent". Several output format can
                be written in the same method call. Defaults to output all codes.
            volSDEF (bool, optional):  Write SDEF definition and tally of solid
                cell for stochastic volume checking. Defaults to False.
            volCARD (bool, optional): Write the CAD calculated volume in the
                cell definition using the VOL card. Defaults to True.
            UCARD (int, optional): Write universe card in the cell definition
                with the specified universe number (if value = 0 Universe card
                is not written). Defaults to None.
            dummyMat (bool, optional): Write dummy material definition card in
               the MCNP output file for all material labels present in the
               model. Dummy material definition is "MX 1001 1". Defaults to False.
            cellCommentFile (bool, optional): Write an additional file with
               comment associated to each CAD cell in the MCNP output file.
               Defaults to False.
            cellSummaryFile (bool, optional): Write an additional file with
               information on the CAD cell translated. Defaults to True.
        """

        write_geometry(
            UniverseBox=self.UniverseBox,
            MetaList=self.MetaList,
            Surfaces=self.Surfaces,
            settings=self.settings,
            options=self.options,
            tolerances=self.tolerances,
            numeric_format=self.numeric_format,
            geometryName=geometryName,
            outFormat=outFormat,
            cellCommentFile=cellCommentFile,
            cellSummaryFile=cellSummaryFile,
            title=title,
            volSDEF=volSDEF,
            volCARD=volCARD,
            UCARD=UCARD,
            dummyMat=dummyMat,
            stepFile=self.stepFile,
        )

        logger.info("End of Monte Carlo code translation phase")

    @classmethod
    def from_json(cls, filename: str):
        """Creates a CadToCsg instance and returns the instance. Populating the
        Options, Tolerance, Settings and NumericFormat attributes from matching
        key names in the JSON. If export_to_csg key is present then this method
        also runs .start() and .export_to_csg() on the instance.

        Args:
            filename str: The filename of the config file. Defaults to "config.json".

        Raises:
            FileNotFoundError: If the config file is not found
            ValueError: If the config JSON file is found to contain an invalid key

        Returns:
            geouned.CadToCsg: returns a geouned CadToCsg class instance.
        """

        if not Path(filename).exists():
            raise FileNotFoundError(f"config file {filename} not found")

        with open(filename) as f:
            config = json.load(f)

        cad_to_csg = cls(stepFile=config["stepFile"])
        for key in config.keys():

            if key in ["stepFile", "export_csg"]:
                pass  # these two keys are used before or after this loop

            elif key == "Tolerances":
                cad_to_csg.tolerances = Tolerances(**config["Tolerances"])

            elif key == "Options":
                cad_to_csg.options = Options(**config["Options"])

            elif key == "NumericFormat":
                cad_to_csg.numeric_format = NumericFormat(**config["NumericFormat"])

            elif key == "Settings":
                cad_to_csg.settings = Settings(**config["Settings"])

            else:
                raise ValueError(
                    f"Invalid key '{key}' found in config file {filename}. Acceptable key names are 'stepFile', 'export_csg', 'Settings', 'Parameters', 'Tolerances' and 'NumericFormat'"
                )

        cad_to_csg.start()
        if "export_csg" in config.keys():
            cad_to_csg.export_csg(**config["export_csg"])
        else:
            cad_to_csg.export_csg()
        return cad_to_csg

    # TODO add tests as set_configuration is not currently tested
    def set_configuration(self, configFile=None):

        if configFile is None:
            return

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(configFile)
        for section in config.sections():
            if section == "Files":
                for key in config["Files"].keys():
                    if key in ("geometryName", "matFile", "title"):
                        self.set(key, config.get("Files", key))

                    elif key == "stepFile":
                        value = config.get("Files", key).strip()
                        lst = value.split()
                        if value[0] in ("(", "[") and value[-1] in ("]", ")"):
                            data = value[1:-1].split(",")
                            data = [x.strip() for x in data]
                            self.set(key, data)
                        elif len(lst) > 1:
                            self.set(key, lst)
                        else:
                            self.set(key, value)

                    elif key == "outFormat":
                        raw = config.get("Files", key).strip()
                        values = tuple(x.strip() for x in raw.split(","))
                        outFormat = []
                        for v in values:
                            if v.lower() == "mcnp":
                                outFormat.append("mcnp")
                            elif v.lower() == "openmc_xml":
                                outFormat.append("openMC_XML")
                            elif v.lower() == "openmc_py":
                                outFormat.append("openMC_PY")
                            elif v.lower() == "serpent":
                                outFormat.append("serpent")
                            elif v.lower() == "phits":
                                outFormat.append("phits")
                        self.set(key, tuple(outFormat))

            elif section == "Parameters":
                for key in config["Parameters"].keys():
                    if key in (
                        "voidGen",
                        "debug",
                        "compSolids",
                        "volSDEF",
                        "volCARD",
                        "dummyMat",
                        "cellSummaryFile",
                        "cellCommentFile",
                        "sort_enclosure",
                    ):
                        self.set(key, config.getboolean("Parameters", key))
                    elif key in (
                        "minVoidSize",
                        "maxSurf",
                        "maxBracket",
                        "startCell",
                        "startSurf",
                    ):
                        self.set(key, config.getint("Parameters", key))
                    elif key in ("exportSolids", "UCARD", "simplify"):
                        self.set(key, config.get("Parameters", key))
                    elif key == "voidMat":
                        value = config.get("Parameters", key).strip()
                        data = value[1:-1].split(",")
                        self.set(key, (int(data[0]), float(data[1]), data[2]))
                    else:
                        value = config.get("Parameters", key).strip()
                        data = value[1:-1].split(",")
                        self.set(key, tuple(map(int, data)))

            elif section == "Options":
                attributes_and_types = get_type_hints(Options())
                for key in config["Options"].keys():
                    if key in attributes_and_types.keys():
                        if attributes_and_types[key] is bool:
                            value = config.getboolean("Options", key)
                        elif (
                            attributes_and_types[key] is float
                            or attributes_and_types[key] is int
                        ):
                            value = config.getfloat("Options", key)
                        setattr(self.options, key, value)

            elif section == "Tolerances":
                attributes_and_types = get_type_hints(Tolerances())
                for key in config["Tolerances"].keys():
                    if key in attributes_and_types.keys():
                        if attributes_and_types[key] is bool:
                            value = config.getboolean("Tolerances", key)
                        elif attributes_and_types[key] is float:
                            value = config.getfloat("Tolerances", key)
                        setattr(self.tolerances, key, value)

            elif section == "MCNP_Numeric_Format":
                attributes_and_types = get_type_hints(NumericFormat())
                PdEntry = False
                for key in config["MCNP_Numeric_Format"].keys():
                    if key in attributes_and_types.keys():
                        value = config.get("MCNP_Numeric_Format", key)
                        setattr(self.numeric_format, key, value)
                        if key == "P_d":
                            PdEntry = True

            else:
                logger.info(f"bad section name : {section}")

        if self.__dict__["geometryName"] == "":
            self.__dict__["geometryName"] = self.__dict__["stepFile"][:-4]

        # TODO see if we can find another way to do this
        if self.options.prnt3PPlane and not PdEntry:
            self.NumericFormat.P_d = "22.15e"

        logger.info(self.__dict__)

    # TODO add tests as set is not currently tested
    def set(self, kwrd, value):

        if kwrd == "stepFile":
            if isinstance(value, (list, tuple)):
                for v in value:
                    if not isinstance(v, str):
                        logger.info(f"elemt in {kwrd} list should be string")
                        return
            elif not isinstance(value, str):
                logger.info(f"{kwrd} should be string or tuple of strings")
                return

        elif kwrd == "UCARD":
            if value == "None":
                value = None
            elif value.isdigit():
                value = int(value)
            else:
                logger.info(f"{kwrd} value should be None or integer")
                return
        elif kwrd == "outFormat":
            if len(value) == 0:
                return
        elif kwrd in ("geometryName", "matFile", "exportSolids"):
            if not isinstance(value, str):
                logger.info(f"{kwrd} value should be str instance")
                return
        elif kwrd in ("cellRange", "voidMat", "voidExclude"):
            if not isinstance(value, (list, tuple)):
                logger.info(f"{kwrd} value should be list or tuple")
                return
        elif kwrd in ("minVoidSize", "maxSurf", "maxBracket", "startCell", "startSurf"):
            if not isinstance(value, int):
                logger.info(f"{kwrd} value should be integer")
                return
        elif kwrd in (
            "voidGen",
            "debug",
            "compSolids",
            "simplifyCTable",
            "volSDEF",
            "volCARD",
            "dummyMat",
            "cellSummaryFile",
            "cellCommentFile",
            "sort_enclosure",
        ):
            if not isinstance(value, bool):
                logger.info(f"{kwrd} value should be boolean")
                return

        self.__dict__[kwrd] = value
        if kwrd == "stepFile" and self.__dict__["geometryName"] == "":
            if isinstance(value, (tuple, list)):
                self.__dict__["geometryName"] == "joined_step_files"
            else:
                self.__dict__["geometryName"] == value[:-4]

    def start(self):

        logger.info("start")
        freecad_version = ".".join(FreeCAD.Version()[:3])
        logger.info(
            f"GEOUNED version {GEOUNED_Version} {GEOUNED_ReleaseDate} \nFreeCAD version {freecad_version}"
        )

        if self.stepFile == "":
            raise ValueError("Cannot run the code. Step file name is missing")

        if isinstance(self.stepFile, (tuple, list)):
            for stp in self.stepFile:
                if not path.isfile(stp):
                    raise FileNotFoundError(f"Step file {stp} not found.\nStop.")
        else:
            if not path.isfile(self.stepFile):
                raise FileNotFoundError(f"Step file {self.stepFile} not found.\nStop.")

        startTime = datetime.now()

        if isinstance(self.stepFile, (list, tuple)):
            step_files = self.stepFile
        else:
            step_files = [self.stepFile]
        MetaChunk = []
        EnclosureChunk = []
        for stp in tqdm(step_files, desc="Loading CAD files"):
            logger.info(f"read step file : {stp}")
            Meta, Enclosure = Load.load_cad(stp, self.settings, self.options)
            MetaChunk.append(Meta)
            EnclosureChunk.append(Enclosure)
        self.MetaList = join_meta_lists(MetaChunk)
        EnclosureList = join_meta_lists(EnclosureChunk)

        logger.info("End of loading phase")
        tempstr1 = str(datetime.now() - startTime)
        logger.info(tempstr1)
        tempTime = datetime.now()

        # Select a specific solid range from original STEP solids
        if self.settings.cellRange:
            self.MetaList = self.MetaList[
                self.settings.cellRange[0] : self.settings.cellRange[1]
            ]

        # export in STEP format solids read from input file
        # terminate excution
        if self.settings.exportSolids != "":
            solids = []
            for m in self.MetaList:
                if m.IsEnclosure:
                    continue
                solids.extend(m.Solids)
            Part.makeCompound(solids).exportStep(self.settings.exportSolids)
            msg = (
                f"Solids exported in file {self.settings.exportSolids}\n"
                "GEOUNED Finish. No solid translation performed."
            )
            raise ValueError(msg)

        # set up Universe
        if EnclosureList:
            self.UniverseBox = get_universe(self.MetaList + EnclosureList)
        else:
            self.UniverseBox = get_universe(self.MetaList)

        self.Surfaces = UF.SurfacesDict(offset=self.settings.startSurf - 1)

        warnSolids = []
        warnEnclosures = []
        coneInfo = dict()
        tempTime0 = datetime.now()
        if not self.options.Facets:

            # decompose all solids in elementary solids (convex ones)
            warningSolidList = decompose_solids(
                self.MetaList,
                self.Surfaces,
                self.UniverseBox,
                self.settings,
                True,
                self.options,
                self.tolerances,
                self.numeric_format,
            )

            # decompose Enclosure solids
            if self.settings.voidGen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    EnclosureList,
                    self.Surfaces,
                    self.UniverseBox,
                    self.settings,
                    False,
                    self.options,
                    self.tolerances,
                    self.numeric_format,
                )

            logger.info("End of decomposition phase")

            # start Building CGS cells phase

            for j, m in enumerate(tqdm(self.MetaList, desc="Translating solid cells")):
                if m.IsEnclosure:
                    continue
                logger.info(f"Building cell: {j+1}")
                cones = Conv.cellDef(
                    m,
                    self.Surfaces,
                    self.UniverseBox,
                    self.options,
                    self.tolerances,
                    self.numeric_format,
                )
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningSolidList:
                    warnSolids.append(m)
                if not m.Solids:
                    logger.info(f"none {j}, {m.__id__}")
                    logger.info(m.Definition)

            if self.options.forceNoOverlap:
                Conv.no_overlapping_cell(self.MetaList, self.Surfaces, self.options)

        else:
            translate(
                self.MetaList,
                self.Surfaces,
                self.UniverseBox,
                self.settings,
                self.options,
                self.tolerances,
            )
            # decompose Enclosure solids
            if self.settings.voidGen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    EnclosureList,
                    self.Surfaces,
                    self.UniverseBox,
                    self.settings,
                    False,
                    self.options,
                    self.tolerances,
                    self.numeric_format,
                )

        tempstr2 = str(datetime.now() - tempTime)
        logger.info(tempstr2)

        #  building enclosure solids

        if self.settings.voidGen and EnclosureList:
            for j, m in enumerate(EnclosureList):
                logger.info(f"Building Enclosure Cell: {j + 1}")
                cones = Conv.cellDef(
                    m,
                    self.Surfaces,
                    self.UniverseBox,
                    self.options,
                    self.tolerances,
                    self.numeric_format,
                )
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningEnclosureList:
                    warnEnclosures.append(m)

        tempTime1 = datetime.now()

        # void generation phase
        MetaVoid = []
        if self.settings.voidGen:
            logger.info("Build Void")
            logger.info(self.settings.voidExclude)
            if not self.settings.voidExclude:
                MetaReduced = self.MetaList
            else:
                MetaReduced = exclude_cells(self.MetaList, self.settings.voidExclude)

            if self.MetaList:
                init = self.MetaList[-1].__id__ - len(EnclosureList)
            else:
                init = 0
            MetaVoid = Void.void_generation(
                MetaReduced,
                EnclosureList,
                self.Surfaces,
                self.UniverseBox,
                self.settings,
                init,
                self.options,
                self.tolerances,
                self.numeric_format,
            )

        # if self.settings.simplify == 'full' and not self.options.forceNoOverlap:
        if self.settings.simplify == "full":
            Surfs = {}
            for lst in self.Surfaces.values():
                for s in lst:
                    Surfs[s.Index] = s

            for c in tqdm(self.MetaList, desc="Simplifying"):
                if c.Definition.level == 0 or c.IsEnclosure:
                    continue
                logger.info(f"simplify cell {c.__id__}")
                Box = UF.get_box(c)
                CT = build_c_table_from_solids(Box, (c.Surfaces, Surfs), option="full")
                c.Definition.simplify(CT)
                c.Definition.clean()
                if type(c.Definition.elements) is bool:
                    logger.info(
                        f"unexpected constant cell {c.__id__} :{c.Definition.elements}"
                    )

        tempTime2 = datetime.now()
        logger.info(f"build Time: {tempTime2} - {tempTime1}")

        logger.info(datetime.now() - startTime)

        cellOffSet = self.settings.startCell - 1
        if EnclosureList and self.settings.sort_enclosure:
            # sort group solid cell / void cell sequence in each for each enclosure
            # if a solid belong to several enclosure, its definition will be written
            # for the highest enclosure level or if same enclosure level in the first
            # enclosure found
            self.MetaList = sort_enclosure(self.MetaList, MetaVoid, cellOffSet)
        else:
            # remove Null Cell and apply cell numbering offset
            deleted = []
            idLabel = {0: 0}
            icount = cellOffSet
            for i, m in enumerate(self.MetaList):
                if m.NullCell or m.IsEnclosure:
                    deleted.append(i)
                    continue

                icount += 1
                m.label = icount
                idLabel[m.__id__] = m.label

            for i in reversed(deleted):
                del self.MetaList[i]

            lineComment = """\
##########################################################
             VOID CELLS
##########################################################"""
            mc = UF.GeounedSolid(None)
            mc.Comments = lineComment
            self.MetaList.append(mc)

            deleted = []
            for i, m in enumerate(MetaVoid):
                if m.NullCell:
                    deleted.append(i)
                    continue
                icount += 1
                m.label = icount
                update_comment(m, idLabel)
            for i in reversed(deleted):
                del MetaVoid[i]

            self.MetaList.extend(MetaVoid)

        print_warning_solids(warnSolids, warnEnclosures)

        # add plane definition to cone
        process_cones(
            self.MetaList,
            coneInfo,
            self.Surfaces,
            self.UniverseBox,
            self.options,
            self.tolerances,
            self.numeric_format,
        )

        logger.info("Process finished")
        logger.info(datetime.now() - startTime)

        logger.info(f"Translation time of solid cells {tempTime1} - {tempTime0}")
        logger.info(f"Translation time of void cells {tempTime2} - {tempTime1}")


def decompose_solids(
    MetaList, Surfaces, UniverseBox, setting, meta, options, tolerances, numeric_format
):
    totsolid = len(MetaList)
    warningSolids = []
    for i, m in enumerate(tqdm(MetaList, desc="Decomposing solids")):
        if meta and m.IsEnclosure:
            continue
        logger.info(f"Decomposing solid: {i + 1}/{totsolid}")
        if setting.debug:
            logger.info(m.Comments)
            if not path.exists("debug"):
                mkdir("debug")
            if m.IsEnclosure:
                m.Solids[0].exportStep(f"debug/origEnclosure_{i}.stp")
            else:
                m.Solids[0].exportStep(f"debug/origSolid_{i}.stp")

        comsolid, err = Decom.SplitSolid(
            Part.makeCompound(m.Solids),
            UniverseBox,
            options,
            tolerances,
            numeric_format,
        )

        if err != 0:
            if not path.exists("Suspicious_solids"):
                mkdir("Suspicious_solids")
            if m.IsEnclosure:
                Part.CompSolid(m.Solids).exportStep(
                    f"Suspicious_solids/Enclosure_original_{i}.stp"
                )
                comsolid.exportStep(f"Suspicious_solids/Enclosure_split_{i}.stp")
            else:
                Part.CompSolid(m.Solids).exportStep(
                    f"Suspicious_solids/Solid_original_{i}.stp"
                )
                comsolid.exportStep(f"Suspicious_solids/Solid_split_{i}.stp")

            warningSolids.append(i)

        if setting.debug:
            if m.IsEnclosure:
                comsolid.exportStep(f"debug/compEnclosure_{i}.stp")
            else:
                comsolid.exportStep(f"debug/compSolid_{i}.stp")
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
        m.set_cad_solid()
        m.update_solids(comsolid.Solids)

    return warningSolids


def update_comment(meta, idLabel):
    if meta.__commentInfo__ is None:
        return
    if meta.__commentInfo__[1] is None:
        return
    newLabel = (idLabel[i] for i in meta.__commentInfo__[1])
    meta.set_comments(Void.void_comment_line((meta.__commentInfo__[0], newLabel)))


def process_cones(
    MetaList, coneInfo, Surfaces, UniverseBox, options, tolerances, numeric_format
):
    cellId = tuple(coneInfo.keys())
    for m in MetaList:
        if m.__id__ not in cellId and not m.Void:
            continue

        if m.Void and m.__commentInfo__ is not None:
            if m.__commentInfo__[1] is None:
                continue
            cones = set()
            for Id in m.__commentInfo__[1]:
                if Id in cellId:
                    cones.update(-x for x in coneInfo[Id])
            Conv.add_cone_plane(
                m.Definition,
                cones,
                Surfaces,
                UniverseBox,
                options,
                tolerances,
                numeric_format,
            )
        elif not m.Void:
            Conv.add_cone_plane(
                m.Definition,
                coneInfo[m.__id__],
                Surfaces,
                UniverseBox,
                options,
                tolerances,
                numeric_format,
            )


def get_universe(MetaList):
    d = 10
    Box = MetaList[0].BoundBox
    xmin = Box.XMin
    xmax = Box.XMax
    ymin = Box.YMin
    ymax = Box.YMax
    zmin = Box.ZMin
    zmax = Box.ZMax
    for m in MetaList[1:]:
        # MIO. This was removed since in HELIAS the enclosure cell is the biggest one
        # if m.IsEnclosure: continue
        xmin = min(m.BoundBox.XMin, xmin)
        xmax = max(m.BoundBox.XMax, xmax)
        ymin = min(m.BoundBox.YMin, ymin)
        ymax = max(m.BoundBox.YMax, ymax)
        zmin = min(m.BoundBox.ZMin, zmin)
        zmax = max(m.BoundBox.ZMax, zmax)

    return FreeCAD.BoundBox(
        FreeCAD.Vector(xmin - d, ymin - d, zmin - d),
        FreeCAD.Vector(xmax + d, ymax + d, zmax + d),
    )


def print_warning_solids(warnSolids, warnEnclosures):

    if warnSolids or warnEnclosures:
        fic = open("Warning_Solids_definition.txt", "w")
    else:
        return

    if warnSolids:
        lines = "Solids :\n"
        for sol in warnSolids:
            lines += "\n"
            lines += f"{sol.label}\n"
            lines += f"{sol.Comments}\n"
            lines += f"{write_mcnp_cell_def(sol.Definition)}\n"
        fic.write(lines)

    if warnEnclosures:
        lines = "Enclosures :\n"
        for sol in warnEnclosures:
            lines += "\n"
            lines += f"{sol.label}\n"
            lines += f"{sol.Comments}\n"
            lines += f"{write_mcnp_cell_def(sol.Definition)}\n"

        fic.write(lines)

    fic.close()


def join_meta_lists(MList):

    newMetaList = MList[0]
    if MList[0]:
        for M in MList[1:]:
            lastID = newMetaList[-1].__id__ + 1
            for i, meta in enumerate(M):
                meta.__id__ = lastID + i
                newMetaList.append(meta)

    return newMetaList


def exclude_cells(MetaList, labelList):
    voidMeta = []
    for m in MetaList:
        if m.IsEnclosure:
            continue
        found = False
        for label in labelList:
            if label in m.Comments:
                found = True
                break
        if not found:
            voidMeta.append(m)

    return voidMeta


def sort_enclosure(MetaList, MetaVoid, offSet=0):

    newList = {}
    for m in MetaVoid:
        if m.EnclosureID in newList.keys():
            newList[m.EnclosureID].append(m)
        else:
            newList[m.EnclosureID] = [m]

    icount = offSet
    idLabel = {0: 0}
    newMeta = []
    for m in MetaList:
        if m.NullCell:
            continue
        if m.IsEnclosure:
            lineComment = f"""##########################################################
             ENCLOSURE {m.EnclosureID}
##########################################################"""
            mc = UF.GeounedSolid(None)
            mc.Comments = lineComment
            newMeta.append(mc)
            for e in newList[m.EnclosureID]:
                if e.NullCell:
                    continue
                icount += 1
                e.label = icount
                idLabel[e.__id__] = e.label
                newMeta.append(e)
            lineComment = f"""##########################################################
            END  ENCLOSURE {m.EnclosureID}
##########################################################"""
            mc = UF.GeounedSolid(None)
            mc.Comments = lineComment
            newMeta.append(mc)

        else:
            icount += 1
            m.label = icount
            idLabel[m.__id__] = m.label
            newMeta.append(m)

    lineComment = """\
##########################################################
             VOID CELLS 
##########################################################"""
    mc = UF.GeounedSolid(None)
    mc.Comments = lineComment
    newMeta.append(mc)

    for v in newList[0]:
        if v.NullCell:
            continue
        icount += 1
        v.label = icount
        idLabel[v.__id__] = v.label
        newMeta.append(v)

    for m in newMeta:
        if not m.Void:
            continue
        if m.IsEnclosure:
            continue
        update_comment(m, idLabel)

    return newMeta
