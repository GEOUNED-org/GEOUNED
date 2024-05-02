# Init file of GEOUNED module
#

# We load the STEP and the materials

import configparser
import typing
from datetime import datetime
from os import mkdir, path
from typing import get_type_hints

import FreeCAD
import Part

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
        stepFile: str,
        title: str = "Geouned conversion",
        geometryName: str = "converted_with_geouned",
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
        options=Options(),
        numeric_format=McnpNumericFormat(),
        tolerances=Tolerances(),
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
        self.options = options
        self.numeric_format = numeric_format
        self.tolerances = tolerances

    @classmethod
    def from_config(cls, filename):

        csg_to_csg = cls()

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(filename)
        for section in config.sections():
            if section == "Files":
                for key in config["Files"].keys():
                    if key in ("geometryName", "matFile", "title"):
                        setattr(csg_to_csg, key, config.get("Files", key))
                    elif key == "stepFile":
                        value = config.get("Files", key).strip()
                        lst = value.split()
                        if value[0] in ("(", "[") and value[-1] in ("]", ")"):
                            data = value[1:-1].split(",")
                            data = [x.strip() for x in data]
                            setattr(csg_to_csg, key, data)
                        elif len(lst) > 1:
                            setattr(csg_to_csg, key, lst)
                        else:
                            setattr(csg_to_csg, key, value)

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
                        setattr(csg_to_csg, key, tuple(outFormat))

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
                        setattr(csg_to_csg, key, config.getboolean("Parameters", key))
                    elif key in (
                        "minVoidSize",
                        "maxSurf",
                        "maxBracket",
                        "startCell",
                        "startSurf",
                    ):
                        setattr(csg_to_csg, key, config.getint("Parameters", key))
                    elif key in ("exportSolids", "UCARD", "simplify"):
                        setattr(csg_to_csg, key, config.get("Parameters", key))
                    elif key == "voidMat":
                        value = config.get("Parameters", key).strip()
                        data = value[1:-1].split(",")
                        setattr(
                            csg_to_csg, key, (int(data[0]), float(data[1]), data[2])
                        )
                    else:
                        value = config.get("Parameters", key).strip()
                        data = value[1:-1].split(",")
                        setattr(csg_to_csg, key, tuple(map(int, data)))

            elif section == "Options":
                option_attribute_names_and_types = get_type_hints(Options)
                key_value_pairs = convert_config_keys_to_dict(
                    config, option_attribute_names_and_types, "Options"
                )
                csg_to_csg.options = Options(**key_value_pairs)

            elif section == "Tolerances":
                option_attribute_names_and_types = get_type_hints(Tolerances)
                key_value_pairs = convert_config_keys_to_dict(
                    config, option_attribute_names_and_types, "Tolerances"
                )
                csg_to_csg.tolerances = Tolerances(**key_value_pairs)

            elif section == "MCNP_Numeric_Format":
                option_attribute_names_and_types = get_type_hints(Tolerances)
                key_value_pairs = convert_config_keys_to_dict(
                    config, option_attribute_names_and_types, "Tolerances"
                )
                # TODO check why special case for p_d this is needed, perhaps it can just be another attribute set by user
                if "P_d" in key_value_pairs.keys():
                    PdEntry = True
                else:
                    PdEntry = False
                csg_to_csg.numeric_format = McnpNumericFormat(**key_value_pairs)

            else:
                print(f"bad section name : {section} found in config file {filename}")

        if csg_to_csg.options.prnt3PPlane and not PdEntry:
            McnpNumericFormat.P_d = "22.15e"

    def start(self):

        print("start")
        FreeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())
        print(
            "GEOUNED version {} {} \nFreeCAD version {}".format(
                GEOUNED_Version, GEOUNED_ReleaseDate, FreeCAD_Version
            )
        )

        if self.__dict__ is None:
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
            MetaChunk = []
            EnclosureChunk = []
            for stp in self.stepFile:
                print(f"read step file : {stp}")
                Meta, Enclosure = Load.load_cad(stp, self.matFile)
                MetaChunk.append(Meta)
                EnclosureChunk.append(Enclosure)
            MetaList = join_meta_lists(MetaChunk)
            EnclosureList = join_meta_lists(EnclosureChunk)
        else:
            print(f"read step file : {self.stepFile}")
            MetaList, EnclosureList = Load.load_cad(
                self.stepFile, self.matFile, self.voidMat, self.compSolids
            )

        print("End of loading phase")
        tempstr1 = str(datetime.now() - startTime)
        print(tempstr1)
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
        if not self.options.Facets:

            # decompose all solids in elementary solids (convex ones)
            warningSolidList = decompose_solids(
                MetaList=MetaList,
                Surfaces=Surfaces,
                UniverseBox=UniverseBox,
                debug=self.debug,
                meta=True,
                nPlaneReverse=self.options.nPlaneReverse,
                splitTolerance=self.options.splitTolerance,
                pln_distance=self.tolerances.pln_distance,
                pln_angle=self.tolerances.pln_angle,
                relativeTol=self.tolerances.relativeTol,
                verbose=self.options.verbose,
                tor_distance=self.tolerances.tor_distance,
                tor_angle=self.tolerances.tor_angle,
                scale_up=self.options.scaleUp,
                cyl_angle=self.tolerances.cyl_angle,
                cyl_distance=self.tolerances.cyl_distance,
            )

            # decompose Enclosure solids
            if self.voidGen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    MetaList=EnclosureList,
                    Surfaces=Surfaces,
                    UniverseBox=UniverseBox,
                    debug=self.debug,
                    meta=False,
                    nPlaneReverse=self.options.nPlaneReverse,
                    splitTolerance=self.options.splitTolerance,
                    pln_distance=self.tolerances.pln_distance,
                    pln_angle=self.tolerances.pln_angle,
                    relativeTol=self.tolerances.relativeTol,
                    verbose=self.options.verbose,
                    tor_distance=self.tolerances.tor_distance,
                    tor_angle=self.tolerances.tor_angle,
                    scale_up=self.options.scaleUp,
                    cyl_angle=self.tolerances.cyl_angle,
                    cyl_distance=self.tolerances.cyl_distance,
                )

            print("End of decomposition phase")

            # start Building CGS cells phase

            for j, m in enumerate(MetaList):
                if m.IsEnclosure:
                    continue
                print("Building cell: ", j + 1)
                cones = Conv.cellDef(
                    m, Surfaces, UniverseBox, self.options.forceCylinder
                )
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningSolidList:
                    warnSolids.append(m)
                if not m.Solids:
                    print("none", j, m.__id__)
                    print(m.Definition)

            if Options.forceNoOverlap:
                Conv.no_overlapping_cell(MetaList, Surfaces)

        else:
            translate(MetaList, Surfaces, UniverseBox, self.debug, self.options.verbose)
            # decompose Enclosure solids
            if self.voidGen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    MataList=EnclosureList,
                    Surfaces=Surfaces,
                    UniverseBox=UniverseBox,
                    debug=self.debug,
                    meta=False,
                    nPlaneReverse=self.options.nPlaneReverse,
                    splitTolerance=self.options.splitTolerance,
                    pln_distance=self.tolerances.pln_distance,
                    pln_angle=self.tolerances.pln_angle,
                    relativeTol=self.tolerances.relativeTol,
                    verbose=self.options.verbose,
                    tor_distance=self.tolerances.tor_distance,
                    tor_angle=self.tolerances.tor_angle,
                    scale_up=self.options.scaleUp,
                    cyl_angle=self.tolerances.cyl_angle,
                    cyl_distance=self.tolerances.cyl_distance,
                )

        tempstr2 = str(datetime.now() - tempTime)
        print(tempstr2)

        #  building enclosure solids

        if self.voidGen and EnclosureList:
            for j, m in enumerate(EnclosureList):
                print("Building Enclosure Cell: ", j + 1)
                cones = Conv.cellDef(
                    m, Surfaces, UniverseBox, self.options.forceCylinder
                )
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningEnclosureList:
                    warnEnclosures.append(m)

        tempTime1 = datetime.now()

        # void generation phase
        MetaVoid = []
        if self.voidGen:
            print("Build Void")
            print(self.voidExclude)
            if not self.voidExclude:
                MetaReduced = MetaList
            else:
                MetaReduced = exclude_cells(MetaList, self.voidExclude)

            if MetaList:
                init = MetaList[-1].__id__ - len(EnclosureList)
            else:
                init = 0
            # TODO perhaps this method should be moved into the CsgToCsg class to avoid passing in so many args
            MetaVoid = Void.void_generation(
                MetaReduced,
                EnclosureList,
                Surfaces,
                UniverseBox,
                self.voidMat,
                self.maxSurf,
                self.maxBracket,
                self.minVoidSize,
                self.simplify,
                self.sort_enclosure,
                init,
                self.options.enlargeBox,
                self.options.verbose,
            )

        # if self.simplify == 'full' and not Options.forceNoOverlap:
        if self.simplify == "full":
            Surfs = {}
            for lst in Surfaces.values():
                for s in lst:
                    Surfs[s.Index] = s

            for c in MetaList:
                if c.Definition.level == 0 or c.IsEnclosure:
                    continue
                print("simplify cell", c.__id__)
                Box = UF.get_box(c)
                CT = build_c_table_from_solids(
                    Box=Box,
                    SurfInfo=(c.Surfaces, Surfs),
                    scale_up=self.options.scale_up,
                    option="full",
                )
                c.Definition.simplify(CT)
                c.Definition.clean()
                if type(c.Definition.elements) is bool:
                    print(
                        f"unexpected constant cell {c.__id__} :{c.Definition.elements}"
                    )

        tempTime2 = datetime.now()
        print("build Time:", tempTime2 - tempTime1)

        print(datetime.now() - startTime)

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

        if self.title == "":
            if isinstance(self.stepFile, str):
                title = self.StepFile
            else:
                title = "; ".join(self.StepFile)
        else:
            title = self.title

        # write outputformat input
        write_geometry(
            UniverseBox=UniverseBox,
            MetaList=MetaList,
            Surfaces=Surfaces,
            title=title,
            volSDEF=self.volSDEF,
            volCARD=self.volCARD,
            UCARD=self.UCARD,
            dummyMat=self.dummyMat,
            geometryName=self.geometryName,
            outFormat=self.outFormat,
            cellCommentFile=self.cellCommentFile,
            cellSummaryFile=self.cellSummaryFile,
            voidGen=self.voidGen,
            matFile=self.matFile,
            voidMat=self.voidMat,
            startCell=self.startCell,
        )

        print("End of MCNP, OpenMC, Serpent and PHITS translation phase")

        print("Process finished")
        print(datetime.now() - startTime)

        print("Translation time of solid cells", tempTime1 - tempTime0)
        print("Translation time of void cells", tempTime2 - tempTime1)


def decompose_solids(
    MetaList,
    Surfaces,
    UniverseBox,
    debug,
    meta,
    nPlaneReverse,
    splitTolerance,
    relativeTol,
    pln_distance,
    pln_angle,
    verbose,
    tor_distance,
    tor_angle,
    scale_up,
    cyl_angle,
    cyl_distance,
):
    totsolid = len(MetaList)
    warningSolids = []
    for i, m in enumerate(MetaList):
        if meta and m.IsEnclosure:
            continue
        print(f"Decomposing solid: {i + 1}/{totsolid} ")
        if debug:
            print(m.Comments)
            if not path.exists("debug"):
                mkdir("debug")
            if m.IsEnclosure:
                m.Solids[0].exportStep(f"debug/origEnclosure_{i}.stp")
            else:
                m.Solids[0].exportStep(f"debug/origSolid_{i}.stp")

        comsolid, err = Decom.split_solid(
            solidShape=Part.makeCompound(m.Solids),
            nPlaneReverse=nPlaneReverse,
            universe_box=UniverseBox,
            splitTolerance=splitTolerance,
            pln_distance=pln_distance,
            pln_angle=pln_angle,
            relativeTol=relativeTol,
            verbose=verbose,
            tor_distance=tor_distance,
            tor_angle=tor_angle,
            scale_up=scale_up,
            cyl_angle=cyl_angle,
            cyl_distance=cyl_distance,
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

        if debug:
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


def convert_config_keys_to_dict(config, option_attribute_names_and_types, key_name):
    config_dict = {}
    for key, _ in config[key_name].items():
        if key in option_attribute_names_and_types.keys():
            if option_attribute_names_and_types[key] is bool:
                value = config.getboolean(key_name)
            elif option_attribute_names_and_types[key] is float:
                value = config.getfloat(key_name)
            elif option_attribute_names_and_types[key] is int:
                value = config.getint(key_name)
            config_dict[key] = value
        else:
            msg = (
                f"{key} was found in the config Options section "
                "but is not an acceptable key name. Acceptable key "
                f"names are {option_attribute_names_and_types}"
            )
            raise ValueError(msg)
    return config_dict
