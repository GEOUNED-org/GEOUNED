import configparser
import json
import logging
import typing
from datetime import datetime
from pathlib import Path
from typing import get_type_hints
from importlib.metadata import version

import FreeCAD
import Part
from tqdm import tqdm

from .code_version import *
from .conversion import cell_definition as Conv
from .cuboid.translate import translate
from .decompose import decom_one as Decom
from .loadfile import load_step as Load
from .utils import functions as UF
from .utils.boolean_solids import build_c_table_from_solids
from .utils.data_classes import NumericFormat, Options, Settings, Tolerances
from .void import void as void
from .write.functions import write_mcnp_cell_def
from .write.write_files import write_geometry

logger = logging.getLogger("general_logger")


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

    @property
    def stepFile(self):
        return self._stepFile

    @stepFile.setter
    def stepFile(self, value: str):
        if not isinstance(value, str):
            raise TypeError(f"geouned.CadToCsg.stepFile should be a str, not a {type(value)}")
        self._stepFile = value

    @property
    def options(self):
        return self._options

    @options.setter
    def options(self, value: Options):
        if not isinstance(value, Options):
            raise TypeError(f"geouned.CadToCsg.options should be an instance of geouned.Options, not a {type(value)}")
        self._options = value

    @property
    def tolerances(self):
        return self._tolerances

    @tolerances.setter
    def tolerances(self, value: tolerances):
        if not isinstance(value, Tolerances):
            raise TypeError(f"geouned.CadToCsg.tolerances should be an instance of geouned.Tolerances, not a {type(value)}")
        self._tolerances = value

    @property
    def numeric_format(self):
        return self._numeric_format

    @numeric_format.setter
    def numeric_format(self, value: numeric_format):
        if not isinstance(value, NumericFormat):
            raise TypeError(
                f"geouned.CadToCsg.numeric_format should be an instance of geouned.NumericFormat, not a {type(value)}"
            )
        self._numeric_format = value

    @property
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self, value: settings):
        if not isinstance(value, Settings):
            raise TypeError(f"geouned.CadToCsg.settings should be an instance of geouned.Settings, not a {type(value)}")
        self._settings = value

    def export_csg(
        self,
        title: str = "Converted with GEOUNED",
        geometryName: str = "csg",
        outFormat: typing.Tuple[str] = (
            "openmc_xml",
            "openmc_py",
            "serpent",
            "phits",
            "mcnp",
        ),
        volSDEF: bool = False,
        volCARD: bool = True,
        UCARD: typing.Union[int, type(None)] = 101,
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
                geometry. Available format are: "mcnp", "openmc_xml",
                "openmc_py", "phits" and "serpent". Several output format can
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

        if not isinstance(UCARD, int) and not isinstance(UCARD, type(None)):
            raise TypeError(f"UCARD should be of type int or None not {type(UCARD)}")

        for arg, arg_str in (
            (volSDEF, "volSDEF"),
            (volCARD, "volCARD"),
            (dummyMat, "dummyMat"),
            (cellCommentFile, "cellCommentFile"),
            (cellSummaryFile, "cellSummaryFile"),
        ):
            if not isinstance(arg, bool):
                raise TypeError(f"{arg} should be of type bool not {type(arg_str)}")

        for arg, arg_str in ((title, "title"), (geometryName, "geometryName")):
            if not isinstance(arg, str):
                raise TypeError(f"{arg} should be of type str not {type(arg_str)}")

        write_geometry(
            UniverseBox=self.UniverseBox,
            MetaList=self.meta_list,
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
                                outFormat.append("openmc_xml")
                            elif v.lower() == "openmc_py":
                                outFormat.append("openmc_py")
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
                        elif attributes_and_types[key] is float or attributes_and_types[key] is int:
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
        logger.info(f"GEOUNED version {version('geouned')} \nFreeCAD version {freecad_version}")

        if self.stepFile == "":
            raise ValueError("Cannot run the code. Step file name is missing")

        if isinstance(self.stepFile, (tuple, list)):
            for stp in self.stepFile:
                if not Path(stp).is_file():
                    raise FileNotFoundError(f"Step file {stp} not found.")
        else:
            if not Path(self.stepFile).is_file():
                raise FileNotFoundError(f"Step file {self.stepFile} not found.")

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
        self.meta_list: typing.List[UF.GeounedSolid] = join_meta_lists(MetaChunk)
        self.enclosure_list: typing.List[UF.GeounedSolid] = join_meta_lists(EnclosureChunk)

        logger.info("End of loading phase")
        tempstr1 = str(datetime.now() - startTime)
        logger.info(tempstr1)
        tempTime = datetime.now()

        # Select a specific solid range from original STEP solids
        if self.settings.cellRange:
            self.meta_list = self.meta_list[self.settings.cellRange[0] : self.settings.cellRange[1]]

        # export in STEP format solids read from input file
        # terminate excution
        if self.settings.exportSolids != "":
            solids = []
            for m in self.meta_list:
                if m.IsEnclosure:
                    continue
                solids.extend(m.Solids)
            Part.makeCompound(solids).exportStep(self.settings.exportSolids)
            msg = f"Solids exported in file {self.settings.exportSolids}\n" "GEOUNED Finish. No solid translation performed."
            raise ValueError(msg)

        # set up Universe
        if self.enclosure_list:
            self.UniverseBox = get_universe(self.meta_list + self.enclosure_list)
        else:
            self.UniverseBox = get_universe(self.meta_list)

        self.Surfaces = UF.SurfacesDict(offset=self.settings.startSurf - 1)

        warnSolids = []
        warnEnclosures = []
        coneInfo = dict()
        tempTime0 = datetime.now()
        if not self.options.Facets:

            # decompose all solids in elementary solids (convex ones)
            warningSolidList = self._decompose_solids(meta=True)

            # decompose Enclosure solids
            if self.settings.voidGen and self.enclosure_list:
                warningEnclosureList = self._decompose_solids(meta=False)

            logger.info("End of decomposition phase")

            # start Building CGS cells phase

            for j, m in enumerate(tqdm(self.meta_list, desc="Translating solid cells")):
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
                Conv.no_overlapping_cell(self.meta_list, self.Surfaces, self.options)

        else:
            translate(
                self.meta_list,
                self.Surfaces,
                self.UniverseBox,
                self.settings,
                self.options,
                self.tolerances,
            )
            # decompose Enclosure solids
            if self.settings.voidGen and self.enclosure_list:
                warningEnclosureList = self._decompose_solids(meta=False)

        tempstr2 = str(datetime.now() - tempTime)
        logger.info(tempstr2)

        #  building enclosure solids

        if self.settings.voidGen and self.enclosure_list:
            for j, m in enumerate(self.enclosure_list):
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
        meta_void = []
        if self.settings.voidGen:
            logger.info("Build Void")
            logger.info(self.settings.voidExclude)
            if not self.settings.voidExclude:
                meta_reduced = self.meta_list
            else:
                meta_reduced = exclude_cells(self.meta_list, self.settings.voidExclude)

            if self.meta_list:
                init = self.meta_list[-1].__id__ - len(self.enclosure_list)
            else:
                init = 0
            meta_void = void.void_generation(
                meta_reduced,
                self.enclosure_list,
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

            for c in tqdm(self.meta_list, desc="Simplifying"):
                if c.Definition.level == 0 or c.IsEnclosure:
                    continue
                logger.info(f"simplify cell {c.__id__}")
                Box = UF.get_box(c)
                CT = build_c_table_from_solids(Box, (c.Surfaces, Surfs), option="full")
                c.Definition.simplify(CT)
                c.Definition.clean()
                if type(c.Definition.elements) is bool:
                    logger.info(f"unexpected constant cell {c.__id__} :{c.Definition.elements}")

        tempTime2 = datetime.now()
        logger.info(f"build Time: {tempTime2} - {tempTime1}")

        logger.info(datetime.now() - startTime)

        cellOffSet = self.settings.startCell - 1
        if self.enclosure_list and self.settings.sort_enclosure:
            # sort group solid cell / void cell sequence in each for each enclosure
            # if a solid belong to several enclosure, its definition will be written
            # for the highest enclosure level or if same enclosure level in the first
            # enclosure found
            self.meta_list = sort_enclosure(self.meta_list, meta_void, cellOffSet)
        else:
            # remove Null Cell and apply cell numbering offset
            deleted = []
            idLabel = {0: 0}
            icount = cellOffSet
            for i, m in enumerate(self.meta_list):
                if m.NullCell or m.IsEnclosure:
                    deleted.append(i)
                    continue

                icount += 1
                m.label = icount
                idLabel[m.__id__] = m.label

            for i in reversed(deleted):
                del self.meta_list[i]

            lineComment = """\
##########################################################
             VOID CELLS
##########################################################"""
            mc = UF.GeounedSolid(None)
            mc.Comments = lineComment
            self.meta_list.append(mc)

            deleted = []
            for i, m in enumerate(meta_void):
                if m.NullCell:
                    deleted.append(i)
                    continue
                icount += 1
                m.label = icount
                update_comment(m, idLabel)
            for i in reversed(deleted):
                del meta_void[i]

            self.meta_list.extend(meta_void)

        print_warning_solids(warnSolids, warnEnclosures)

        # add plane definition to cone
        process_cones(
            self.meta_list,
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

    def _decompose_solids(self, meta: bool):

        if meta:
            meta_list = self.meta_list
            description = "Decomposing solids"
        else:
            meta_list = self.enclosure_list
            description = "Decomposing enclosure solids"

        totsolid = len(meta_list)
        warningSolids = []
        for i, m in enumerate(tqdm(meta_list, desc=description)):
            if meta and m.IsEnclosure:
                continue
            logger.info(f"Decomposing solid: {i + 1}/{totsolid}")
            if self.settings.debug:
                logger.info(m.Comments)
                Path('debug').mkdir(parents=True, exist_ok=True)
                if m.IsEnclosure:
                    m.Solids[0].exportStep(f"debug/origEnclosure_{i}.stp")
                else:
                    m.Solids[0].exportStep(f"debug/origSolid_{i}.stp")

            comsolid, err = Decom.SplitSolid(
                Part.makeCompound(m.Solids),
                self.UniverseBox,
                self.options,
                self.tolerances,
                self.numeric_format,
            )

            if err != 0:
                Path('suspicious_solids').mkdir(parents=True, exist_ok=True)
                if m.IsEnclosure:
                    Part.CompSolid(m.Solids).exportStep(f"Suspicious_solids/Enclosure_original_{i}.stp")
                    comsolid.exportStep(f"Suspicious_solids/Enclosure_split_{i}.stp")
                else:
                    Part.CompSolid(m.Solids).exportStep(f"Suspicious_solids/Solid_original_{i}.stp")
                    comsolid.exportStep(f"Suspicious_solids/Solid_split_{i}.stp")

                warningSolids.append(i)

            if self.settings.debug:
                if m.IsEnclosure:
                    comsolid.exportStep(f"debug/compEnclosure_{i}.stp")
                else:
                    comsolid.exportStep(f"debug/compSolid_{i}.stp")
            self.Surfaces.extend(
                Decom.extract_surfaces(
                    comsolid,
                    "All",
                    self.UniverseBox,
                    self.options,
                    self.tolerances,
                    self.numeric_format,
                    MakeObj=True,
                ),
                self.options,
                self.tolerances,
                self.numeric_format,
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
    meta.set_comments(void.void_comment_line((meta.__commentInfo__[0], newLabel)))


def process_cones(MetaList, coneInfo, Surfaces, UniverseBox, options, tolerances, numeric_format):
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

    solids_logger = logging.getLogger("solids_logger")

    if warnSolids or warnEnclosures:
        pass
    else:
        return

    if warnSolids:
        lines = "Solids :\n"
        for sol in warnSolids:
            lines += "\n"
            lines += f"{sol.label}\n"
            lines += f"{sol.Comments}\n"
            lines += f"{write_mcnp_cell_def(sol.Definition)}\n"
        solids_logger.info(lines)

    if warnEnclosures:
        lines = "Enclosures :\n"
        for sol in warnEnclosures:
            lines += "\n"
            lines += f"{sol.label}\n"
            lines += f"{sol.Comments}\n"
            lines += f"{write_mcnp_cell_def(sol.Definition)}\n"

        solids_logger.info(lines)


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


def sort_enclosure(MetaList, meta_void, offSet=0):

    newList = {}
    for m in meta_void:
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
