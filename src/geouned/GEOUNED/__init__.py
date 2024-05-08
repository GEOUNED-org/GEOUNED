# Init file of GEOUNED module
#

# We load the STEP and the materials

import configparser
import logging
import typing
from datetime import datetime
from os import mkdir, path

import FreeCAD
import Part
from tqdm import tqdm

from .CodeVersion import *
from .Conversion import CellDefinition as Conv
from .Cuboid.translate import translate
from .Decompose import Decom_one as Decom
from .LoadFile import LoadSTEP as Load
from .Utils import Functions as UF
from .Utils.BooleanSolids import build_c_table_from_solids
from .Utils.Options.Classes import McnpNumericFormat, Options, Tolerances
from .Void import Void as Void
from .Write.Functions import write_mcnp_cell_def
from .Write.WriteFiles import write_geometry

logger = logging.getLogger(__name__)


class CadToCsg:
    """Base class for the conversion of CAD to CSG models

    Args:
        title (str, optional): Title of the model. Defaults to "Geouned
            conversion".
        stepFile (str, optional): Name of the CAD file (in STEP format) to
            be converted. Defaults to "".
        geometryName (str, optional): Base name of the output file(s).
            Defaults to "".
        matFile (str, optional): _description_. Defaults to "".
        outFormat (typing.Tuple[str], optional): Format for the output
            geometry. Available format are: mcnp, openMC_XML, openMC_PY,
            phits and serpent. Several output format can be written in the
            same geouned run. Defaults to ("mcnp",).
        voidGen (bool, optional): Generate voids of the geometry. Defaults
            to True.
        debug (bool, optional): Write step files of original and decomposed
            solids, for each solid in the STEP file. Defaults to False.
        compSolids (bool, optional): Join subsolids of STEP file as a single
            compound solid. Step files generated with SpaceClaim have not
            exactly the same level of solids as FreeCAD. It may a happened
            that solids defined has separated solids are read by FreeCAD
            as a single compound solid (and will produce only one MCNP
            cell). In this case compSolids should be set to False. Defaults
            to True.
        volSDEF (bool, optional): Write SDEF definition and tally of solid
            cell for stochastic volume checking. Defaults to False.
        dummyMat (bool, optional): Write dummy material definition card in
            the MCNP output file for all material labels present in the
            model. Dummy material definition is "MX 1001 1". Defaults to
            False.
        volCARD (bool, optional): Write the CAD calculated volume in the
            cell definition using the VOL card. Defaults to True.
        UCARD (_type_, optional): Write universe card in the cell definition
            with the specified universe number (if value = 0 Universe card
            is not written). Defaults to None.
        simplify (str, optional): Simplify the cell definition considering
            relative surfaces position and using Boolean logics. Available
            options are: "no" no optimization, "void" only void cells are
            simplified. Algorithm is faster but the simplification is not
            optimal. "voidfull" : only void cells are simplified with the
            most optimal algorithm. The time of the conversion can be
            multiplied by 5 or more. "full" : all the cells (solids and
            voids) are simplified. Defaults to "No".
        cellRange (list, optional): Range of cell to be converted (only one
            range is allowed, e.g [100,220]). Default all solids are
            converted. Defaults to [].
        exportSolids (str, optional): Export CAD solid after reading.
            The execution is stopped after export, the translation is not
            carried out. Defaults to "".
        minVoidSize (float, optional): Minimum size of the edges of the
            void cell. Units are in mm. Defaults to 200.0.
        maxBracket (int, optional): Maximum number of brackets (solid
            complementary) allowed in void cell definition. Defaults to 30.
        voidMat (list, optional): Assign a material defined by the user
            instead of void for cells without material definition and the
            cells generated in the automatic void generation. The format
            is a 3 valued tuple (mat_label, mat_density, mat_description).
            Example (100,1e-3,'Air assigned to Void'). Defaults to [].
        voidExclude (list, optional): #TODO see issue 87. Defaults to [].
        startCell (int, optional): Starting cell numbering label. Defaults to 1.
        startSurf (int, optional): Starting surface numbering label. Defaults to 1.
        cellCommentFile (bool, optional): Write an additional file with
            comment associated to each CAD cell in the MCNP output file.
            Defaults to False.
        cellSummaryFile (bool, optional): Write an additional file with
            information on the CAD cell translated. Defaults to True.
        sort_enclosure (bool, optional): If enclosures are defined in the
            CAD models, the voids cells of the enclosure will be located in
            the output file in the same location where the enclosure solid
            is located in the CAD solid tree.. Defaults to False.
    """

    def __init__(
        self,
        title: str = "Geouned conversion",
        stepFile: str = "",
        geometryName: str = "",
        matFile: str = "",
        outFormat: typing.Tuple[str] = ("mcnp",),
        voidGen: bool = True,
        debug: bool = False,
        compSolids: bool = True,
        volSDEF: bool = False,
        dummyMat: bool = False,
        volCARD: bool = True,
        UCARD=None,
        simplify: str = "No",
        cellRange=[],
        exportSolids: str = "",
        minVoidSize: float = 200.0,  # units mm
        maxSurf: int = 50,
        maxBracket: int = 30,
        voidMat=[],
        voidExclude=[],
        startCell: int = 1,
        startSurf: int = 1,
        cellCommentFile: bool = False,
        cellSummaryFile: bool = True,
        sort_enclosure: bool = False,
    ):

        self.title = title
        self.stepFile = stepFile
        self.geometryName = geometryName
        self.matFile = matFile
        self.outFormat = outFormat
        self.voidGen = voidGen
        self.debug = debug
        self.compSolids = compSolids
        self.volSDEF = volSDEF
        self.dummyMat = dummyMat
        self.volCARD = volCARD
        self.UCARD = UCARD
        self.simplify = simplify
        self.cellRange = cellRange
        self.exportSolids = exportSolids
        self.minVoidSize = minVoidSize
        self.maxSurf = maxSurf
        self.maxBracket = maxBracket
        self.voidMat = voidMat
        self.voidExclude = voidExclude
        self.startCell = startCell
        self.startSurf = startSurf
        self.cellCommentFile = cellCommentFile
        self.cellSummaryFile = cellSummaryFile
        self.sort_enclosure = sort_enclosure

        Options.set_default_attribute()
        McnpNumericFormat.set_default_attribute()
        Tolerances.set_default_attribute()

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
                for key in config["Options"].keys():
                    if key in Options.default_values.keys():
                        if Options.type_dict[key] is bool:
                            Options.set_attribute(
                                key, config.getboolean("Options", key)
                            )
                        elif (
                            Options.type_dict[key] is float
                            or Options.type_dict[key] is int
                        ):
                            Options.set_attribute(key, config.getfloat("Options", key))

            elif section == "Tolerances":
                for key in config["Tolerances"].keys():
                    eqvKey = Tolerances.KwrdEquiv[key]
                    if eqvKey in Tolerances.default_values.keys():
                        if Tolerances.type_dict[eqvKey] is bool:
                            Tolerances.set_attribute(
                                eqvKey, config.getboolean("Tolerances", key)
                            )
                        elif Tolerances.type_dict[eqvKey] is float:
                            Tolerances.set_attribute(
                                eqvKey, config.getfloat("Tolerances", key)
                            )

            elif section == "MCNP_Numeric_Format":
                PdEntry = False
                for key in config["MCNP_Numeric_Format"].keys():
                    if key in McnpNumericFormat.default_values.keys():
                        McnpNumericFormat.set_attribute(
                            key, config.get("MCNP_Numeric_Format", key)
                        )
                        if key == "P_d":
                            PdEntry = True

            else:
                logger.info(f"bad section name : {section}")

        if self.__dict__["geometryName"] == "":
            self.__dict__["geometryName"] = self.__dict__["stepFile"][:-4]

        if Options.prnt3PPlane and not PdEntry:
            McnpNumericFormat.P_d = "22.15e"

        logger.info(self.__dict__)

    def set(self, kwrd, value):

        if kwrd in McnpNumericFormat.default_values.keys():
            McnpNumericFormat.set_attribute(kwrd, value)
            return
        elif kwrd in Tolerances.default_values.keys():
            Tolerances.set_attribute(kwrd, value)
            return
        elif kwrd in Options.default_values.keys():
            Options.set_attribute(kwrd, value)
            return
        elif kwrd not in self.__dict__.keys():
            logger.info(f"Bad entry : {kwrd}")
            return

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
        FreeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())
        logger.info(
            "GEOUNED version {} {} \nFreeCAD version {}".format(
                GEOUNED_Version, GEOUNED_ReleaseDate, FreeCAD_Version
            )
        )

        code_setting = self.__dict__
        if code_setting is None:
            raise ValueError("Cannot run the code. Input are missing")
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
        for stp in tqdm(step_files, desc="loading CAD files"):
            logger.info(f"read step file : {stp}")
            Meta, Enclosure = Load.load_cad(stp, self.matFile)
            MetaChunk.append(Meta)
            EnclosureChunk.append(Enclosure)
        MetaList = join_meta_lists(MetaChunk)
        EnclosureList = join_meta_lists(EnclosureChunk)

        logger.info("End of loading phase")
        tempstr1 = str(datetime.now() - startTime)
        logger.info(tempstr1)
        tempTime = datetime.now()

        # Select a specific solid range from original STEP solids
        if self.cellRange:
            MetaList = MetaList[self.cellRange[0] : self.cellRange[1]]

        # export in STEP format solids read from input file
        # terminate excution
        if self.exportSolids != "":
            solids = []
            for m in MetaList:
                if m.IsEnclosure:
                    continue
                solids.extend(m.Solids)
            Part.makeCompound(solids).exportStep(self.exportSolids)
            msg = (
                f"Solids exported in file {self.exportSolids}\n"
                "GEOUNED Finish. No solid translation performed."
            )
            raise ValueError(msg)

        # set up Universe
        if EnclosureList:
            UniverseBox = get_universe(MetaList + EnclosureList)
        else:
            UniverseBox = get_universe(MetaList)

        Surfaces = UF.SurfacesDict(offset=self.startSurf - 1)

        warnSolids = []
        warnEnclosures = []
        coneInfo = dict()
        tempTime0 = datetime.now()
        if not Options.Facets:

            # decompose all solids in elementary solids (convex ones)
            warningSolidList = decompose_solids(
                MetaList, Surfaces, UniverseBox, code_setting, True
            )

            # decompose Enclosure solids
            if self.voidGen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    EnclosureList, Surfaces, UniverseBox, code_setting, False
                )

            logger.info("End of decomposition phase")

            # start Building CGS cells phase

            for j, m in enumerate(MetaList):
                if m.IsEnclosure:
                    continue
                logger.info(f"Building cell: {j+1}")
                cones = Conv.cellDef(m, Surfaces, UniverseBox)
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningSolidList:
                    warnSolids.append(m)
                if not m.Solids:
                    logger.info(f"none {j}, {m.__id__}")
                    logger.info(m.Definition)

            if Options.forceNoOverlap:
                Conv.no_overlapping_cell(MetaList, Surfaces)

        else:
            translate(MetaList, Surfaces, UniverseBox, code_setting)
            # decompose Enclosure solids
            if self.voidGen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    EnclosureList, Surfaces, UniverseBox, code_setting, False
                )

        tempstr2 = str(datetime.now() - tempTime)
        logger.info(tempstr2)

        #  building enclosure solids

        if self.voidGen and EnclosureList:
            for j, m in enumerate(EnclosureList):
                logger.info(f"Building Enclosure Cell: {j + 1}")
                cones = Conv.cellDef(m, Surfaces, UniverseBox)
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningEnclosureList:
                    warnEnclosures.append(m)

        tempTime1 = datetime.now()

        # void generation phase
        MetaVoid = []
        if self.voidGen:
            logger.info("Build Void")
            logger.info(self.voidExclude)
            if not self.voidExclude:
                MetaReduced = MetaList
            else:
                MetaReduced = exclude_cells(MetaList, self.voidExclude)

            if MetaList:
                init = MetaList[-1].__id__ - len(EnclosureList)
            else:
                init = 0
            MetaVoid = Void.void_generation(
                MetaReduced, EnclosureList, Surfaces, UniverseBox, code_setting, init
            )

        # if code_setting['simplify'] == 'full' and not Options.forceNoOverlap:
        if self.simplify == "full":
            Surfs = {}
            for lst in Surfaces.values():
                for s in lst:
                    Surfs[s.Index] = s

            for c in tqdm(MetaList, desc="simplifying"):
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

        cellOffSet = self.startCell - 1
        if EnclosureList and self.sort_enclosure:
            # sort group solid cell / void cell sequence in each for each enclosure
            # if a solid belong to several enclosure, its definition will be written
            # for the highest enclosure level or if same enclosure level in the first
            # enclosure found
            MetaList = sort_enclosure(MetaList, MetaVoid, cellOffSet)
        else:
            # remove Null Cell and apply cell numbering offset
            deleted = []
            idLabel = {0: 0}
            icount = cellOffSet
            for i, m in enumerate(MetaList):
                if m.NullCell or m.IsEnclosure:
                    deleted.append(i)
                    continue

                icount += 1
                m.label = icount
                idLabel[m.__id__] = m.label

            for i in reversed(deleted):
                del MetaList[i]

            lineComment = """\
##########################################################
             VOID CELLS
##########################################################"""
            mc = UF.GeounedSolid(None)
            mc.Comments = lineComment
            MetaList.append(mc)

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

            MetaList.extend(MetaVoid)

        print_warning_solids(warnSolids, warnEnclosures)

        # add plane definition to cone
        process_cones(MetaList, coneInfo, Surfaces, UniverseBox)

        # write outputformat input
        write_geometry(UniverseBox, MetaList, Surfaces, code_setting)

        logger.info("End of MCNP, OpenMC, Serpent and PHITS translation phase")

        logger.info("Process finished")
        logger.info(datetime.now() - startTime)

        logger.info(f"Translation time of solid cells {tempTime1} - {tempTime0}")
        logger.info(f"Translation time of void cells {tempTime2} - {tempTime1}")


def decompose_solids(MetaList, Surfaces, UniverseBox, setting, meta):
    totsolid = len(MetaList)
    warningSolids = []
    for i, m in enumerate(tqdm(MetaList, desc="decomposing solids")):
        if meta and m.IsEnclosure:
            continue
        logger.info(f"Decomposing solid: {i + 1}/{totsolid} ")
        if setting["debug"]:
            logger.info(m.Comments)
            if not path.exists("debug"):
                mkdir("debug")
            if m.IsEnclosure:
                m.Solids[0].exportStep(f"debug/origEnclosure_{i}.stp")
            else:
                m.Solids[0].exportStep(f"debug/origSolid_{i}.stp")

        comsolid, err = Decom.SplitSolid(Part.makeCompound(m.Solids), UniverseBox)

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

        if setting["debug"]:
            if m.IsEnclosure:
                comsolid.exportStep(f"debug/compEnclosure_{i}.stp")
            else:
                comsolid.exportStep(f"debug/compSolid_{i}.stp")
        Surfaces.extend(
            Decom.extract_surfaces(comsolid, "All", UniverseBox, MakeObj=True)
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


def process_cones(MetaList, coneInfo, Surfaces, UniverseBox):
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
            Conv.add_cone_plane(m.Definition, cones, Surfaces, UniverseBox)
        elif not m.Void:
            Conv.add_cone_plane(m.Definition, coneInfo[m.__id__], Surfaces, UniverseBox)


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
