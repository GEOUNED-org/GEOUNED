# Init file of GEOUNED module
#

# We load the STEP and the materials

import configparser
import typing
from datetime import datetime
from os import mkdir, path

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
        title: str = "Geouned conversion",
        step_file: str = "",
        geometry_name: str = "",
        mat_file: str = "",
        out_format: typing.Tuple[str] = ("mcnp",),
        void_gen: bool = True,
        debug: bool = False,
        comp_solids: bool = True,
        vol_sdef: bool = False,
        dummy_mat: bool = False,
        vol_card: bool = True,
        u_card=None,
        simplify: str = "No",
        cell_range=[],
        export_solids: str = "",
        min_void_size: float = 200.0,  # units mm
        max_surf: int = 50,
        max_bracket: int = 30,
        void_mat=[],
        void_exclude=[],
        start_cell: int = 1,
        start_surf: int = 1,
        cell_comment_file: bool = False,
        cell_summary_file: bool = True,
        sort_enclosure: bool = False,
    ):

        self.title = title
        self.step_file = step_file
        self.geometry_name = geometry_name
        self.mat_file = mat_file
        self.out_format = out_format
        self.void_gen = void_gen
        self.debug = debug
        self.comp_solids = comp_solids
        self.vol_sdef = vol_sdef
        self.dummy_mat = dummy_mat
        self.vol_card = vol_card
        self.u_card = u_card
        self.simplify = simplify
        self.cell_range = cell_range
        self.export_solids = export_solids
        self.min_void_size = min_void_size
        self.max_surf = max_surf
        self.max_bracket = max_bracket
        self.void_mat = void_mat
        self.void_exclude = void_exclude
        self.start_cell = start_cell
        self.start_surf = start_surf
        self.cell_comment_file = cell_comment_file
        self.cell_summary_file = cell_summary_file
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
                    if key in ("geometry_name", "mat_file", "title"):
                        self.set(key, config.get("Files", key))

                    elif key == "step_file":
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

                    elif key == "out_format":
                        raw = config.get("Files", key).strip()
                        values = tuple(x.strip() for x in raw.split(","))
                        out_format = []
                        for v in values:
                            if v.lower() == "mcnp":
                                out_format.append("mcnp")
                            elif v.lower() == "openmc_xml":
                                out_format.append("openMC_XML")
                            elif v.lower() == "openmc_py":
                                out_format.append("openMC_PY")
                            elif v.lower() == "serpent":
                                out_format.append("serpent")
                            elif v.lower() == "phits":
                                out_format.append("phits")
                        self.set(key, tuple(out_format))

            elif section == "Parameters":
                for key in config["Parameters"].keys():
                    if key in (
                        "void_gen",
                        "debug",
                        "comp_solids",
                        "vol_sdef",
                        "vol_card",
                        "dummy_mat",
                        "cell_summary_file",
                        "cell_comment_file",
                        "sort_enclosure",
                    ):
                        self.set(key, config.getboolean("Parameters", key))
                    elif key in (
                        "min_void_size",
                        "max_surf",
                        "max_bracket",
                        "start_cell",
                        "start_surf",
                    ):
                        self.set(key, config.getint("Parameters", key))
                    elif key in ("export_solids", "u_card", "simplify"):
                        self.set(key, config.get("Parameters", key))
                    elif key == "void_mat":
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
                print(f"bad section name : {section}")

        if self.__dict__["geometry_name"] == "":
            self.__dict__["geometry_name"] = self.__dict__["step_file"][:-4]

        if Options.prnt3PPlane and not PdEntry:
            McnpNumericFormat.P_d = "22.15e"

        print(self.__dict__)

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
            print(f"Bad entry : {kwrd}")
            return

        if kwrd == "step_file":
            if isinstance(value, (list, tuple)):
                for v in value:
                    if not isinstance(v, str):
                        print(f"elemt in {kwrd} list should be string")
                        return
            elif not isinstance(value, str):
                print(f"{kwrd} should be string or tuple of strings")
                return

        elif kwrd == "u_card":
            if value == "None":
                value = None
            elif value.isdigit():
                value = int(value)
            else:
                print(f"{kwrd} value should be None or integer")
                return
        elif kwrd == "out_format":
            if len(value) == 0:
                return
        elif kwrd in ("geometry_name", "mat_file", "export_solids"):
            if not isinstance(value, str):
                print(f"{kwrd} value should be str instance")
                return
        elif kwrd in ("cell_range", "void_mat", "void_exclude"):
            if not isinstance(value, (list, tuple)):
                print(f"{kwrd} value should be list or tuple")
                return
        elif kwrd in (
            "min_void_size",
            "max_surf",
            "max_bracket",
            "start_cell",
            "start_surf",
        ):
            if not isinstance(value, int):
                print(f"{kwrd} value should be integer")
                return
        elif kwrd in (
            "void_gen",
            "debug",
            "comp_solids",
            "simplifyCTable",
            "vol_sdef",
            "vol_card",
            "dummy_mat",
            "cell_summary_file",
            "cell_comment_file",
            "sort_enclosure",
        ):
            if not isinstance(value, bool):
                print(f"{kwrd} value should be boolean")
                return

        self.__dict__[kwrd] = value
        if kwrd == "step_file" and self.__dict__["geometry_name"] == "":
            if isinstance(value, (tuple, list)):
                self.__dict__["geometry_name"] == "joined_step_files"
            else:
                self.__dict__["geometry_name"] == value[:-4]

    def start(self):

        print("start")
        FreeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())
        print(
            "GEOUNED version {} {} \nFreeCAD version {}".format(
                GEOUNED_Version, GEOUNED_ReleaseDate, FreeCAD_Version
            )
        )

        code_setting = self.__dict__
        if code_setting is None:
            raise ValueError("Cannot run the code. Input are missing")
        if self.step_file == "":
            raise ValueError("Cannot run the code. Step file name is missing")

        if isinstance(self.step_file, (tuple, list)):
            for stp in self.step_file:
                if not path.isfile(stp):
                    raise FileNotFoundError(f"Step file {stp} not found.\nStop.")
        else:
            if not path.isfile(self.step_file):
                raise FileNotFoundError(f"Step file {self.step_file} not found.\nStop.")

        startTime = datetime.now()

        if isinstance(self.step_file, (list, tuple)):
            MetaChunk = []
            EnclosureChunk = []
            for stp in self.step_file:
                print(f"read step file : {stp}")
                Meta, Enclosure = Load.load_cad(stp, self.mat_file)
                MetaChunk.append(Meta)
                EnclosureChunk.append(Enclosure)
            MetaList = join_meta_lists(MetaChunk)
            EnclosureList = join_meta_lists(EnclosureChunk)
        else:
            print(f"read step file : {self.step_file}")
            MetaList, EnclosureList = Load.load_cad(
                self.step_file, self.mat_file, self.void_mat, self.comp_solids
            )

        print("End of loading phase")
        tempstr1 = str(datetime.now() - startTime)
        print(tempstr1)
        tempTime = datetime.now()

        # Select a specific solid range from original STEP solids
        if self.cell_range:
            MetaList = MetaList[self.cell_range[0] : self.cell_range[1]]

        # export in STEP format solids read from input file
        # terminate excution
        if self.export_solids != "":
            solids = []
            for m in MetaList:
                if m.IsEnclosure:
                    continue
                solids.extend(m.Solids)
            Part.makeCompound(solids).exportStep(self.export_solids)
            msg = (
                f"Solids exported in file {self.export_solids}\n"
                "GEOUNED Finish. No solid translation performed."
            )
            raise ValueError(msg)

        # set up Universe
        if EnclosureList:
            UniverseBox = get_universe(MetaList + EnclosureList)
        else:
            UniverseBox = get_universe(MetaList)

        Surfaces = UF.SurfacesDict(offset=self.start_surf - 1)

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
            if self.void_gen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    EnclosureList, Surfaces, UniverseBox, code_setting, False
                )

            print("End of decomposition phase")

            # start Building CGS cells phase

            for j, m in enumerate(MetaList):
                if m.IsEnclosure:
                    continue
                print("Building cell: ", j + 1)
                cones = Conv.cellDef(m, Surfaces, UniverseBox)
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
            translate(MetaList, Surfaces, UniverseBox, code_setting)
            # decompose Enclosure solids
            if self.void_gen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    EnclosureList, Surfaces, UniverseBox, code_setting, False
                )

        tempstr2 = str(datetime.now() - tempTime)
        print(tempstr2)

        #  building enclosure solids

        if self.void_gen and EnclosureList:
            for j, m in enumerate(EnclosureList):
                print("Building Enclosure Cell: ", j + 1)
                cones = Conv.cellDef(m, Surfaces, UniverseBox)
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningEnclosureList:
                    warnEnclosures.append(m)

        tempTime1 = datetime.now()

        # void generation phase
        MetaVoid = []
        if self.void_gen:
            print("Build Void")
            print(self.void_exclude)
            if not self.void_exclude:
                MetaReduced = MetaList
            else:
                MetaReduced = exclude_cells(MetaList, self.void_exclude)

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

            for c in MetaList:
                if c.Definition.level == 0 or c.IsEnclosure:
                    continue
                print("simplify cell", c.__id__)
                Box = UF.get_box(c)
                CT = build_c_table_from_solids(Box, (c.Surfaces, Surfs), option="full")
                c.Definition.simplify(CT)
                c.Definition.clean()
                if type(c.Definition.elements) is bool:
                    print(
                        f"unexpected constant cell {c.__id__} :{c.Definition.elements}"
                    )

        tempTime2 = datetime.now()
        print("build Time:", tempTime2 - tempTime1)

        print(datetime.now() - startTime)

        cellOffSet = self.start_cell - 1
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

        print("End of MCNP, OpenMC, Serpent and PHITS translation phase")

        print("Process finished")
        print(datetime.now() - startTime)

        print("Translation time of solid cells", tempTime1 - tempTime0)
        print("Translation time of void cells", tempTime2 - tempTime1)


def decompose_solids(MetaList, Surfaces, UniverseBox, setting, meta):
    totsolid = len(MetaList)
    warningSolids = []
    for i, m in enumerate(MetaList):
        if meta and m.IsEnclosure:
            continue
        print(f"Decomposing solid: {i + 1}/{totsolid} ")
        if setting["debug"]:
            print(m.Comments)
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
