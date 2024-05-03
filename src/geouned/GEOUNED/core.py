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

from .CodeVersion import GEOUNED_Version, GEOUNED_ReleaseDate
from .Conversion import CellDefinition as Conv
from .Cuboid.translate import translate
from .Decompose import Decom_one as Decom
from .LoadFile import LoadSTEP as Load
from .Utils import Functions as UF
from .Utils.BooleanSolids import build_c_table_from_solids
from .Utils.Classes import NumericFormat, Options, Tolerances, Settings
from .Void import Void as Void
from .Write.Functions import write_mcnp_cell_def
from .Write.WriteFiles import write_geometry


class CadToCsg:
    """Base class for the conversion of CAD to CSG models"""

    def __init__(
        self,
        stepFile: str,
        settings=Settings(),
        options=Options(),
        numeric_format=NumericFormat(),
        tolerances=Tolerances(),
    ):
        self.stepFile = stepFile
        self.settings = settings
        self.options = options
        self.numeric_format = numeric_format
        self.tolerances = tolerances

        # set internally by the class
        self.UniverseBox = None
        self.Surfaces = None
        self.MetaList = None

    @classmethod
    def from_config(cls, filename: str='config.ini') -> 'geouned.CadToCsg':

        csg_to_csg = cls()

        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(filename)
        for section in config.sections():
            if section == "Files":
                for key in config["Files"].keys():
                    if key in ("matFile"):
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
                                outFormat.append("openmc_xml")
                            elif v.lower() == "openmc_py":
                                outFormat.append("openmc_py")
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
                        "sort_enclosure",
                    ):
                        setattr(csg_to_csg, key, config.getboolean("Parameters", key))
                    elif key in (
                        "minVoidSize",
                        "max_surf",
                        "max_bracket",
                        "startCell",
                        "startSurf",
                    ):
                        setattr(csg_to_csg, key, config.getint("Parameters", key))
                    elif key in ("exportSolids", "simplify"):
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
                csg_to_csg.numeric_format = NumericFormat(**key_value_pairs)

            else:
                print(f"bad section name : {section} found in config file {filename}")

        if csg_to_csg.options.prnt3PPlane and not PdEntry:
            csg_to_csg.P_d = "22.15e"

    def start(self):

        if self.options.verbose:
            print("started converting")
            FreeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())
            print(f"GEOUNED version {GEOUNED_Version} {GEOUNED_ReleaseDate} \nFreeCAD version {FreeCAD_Version}")

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
                print(f"reading step file : {stp}")
                Meta, Enclosure = Load.load_cad(
                    filename=stp,
                    mat_filename=self.matFile,
                    default_mat=[],
                    delLastNumber=self.options.delLastNumber,
                    comp_solids=True
                    )
                MetaChunk.append(Meta)
                EnclosureChunk.append(Enclosure)
            self.MetaList = join_meta_lists(MetaChunk)
            EnclosureList = join_meta_lists(EnclosureChunk)
        else:
            print(f"reading step file : {self.stepFile}")
            self.MetaList, EnclosureList = Load.load_cad(
                filename=self.stepFile, mat_filename=self.settings.matFile,
                default_mat=self.settings.voidMat, delLastNumber=self.options.delLastNumber,
                comp_solids=self.settings.compSolids
            )

        if self.options.verbose:
            tempstr1 = str(datetime.now() - startTime)
            print(f"Loading STEP file(s)phase, time = {tempstr1}")
            tempTime = datetime.now()

        # Select a specific solid range from original STEP solids
        if self.settings.cellRange:
            self.MetaList = self.MetaList[self.settings.cellRange[0] : self.settings.cellRange[1]]

        # export in STEP format solids read from input file
        # terminate excution
        if self.settings.exportSolids != "":
            solids = []
            for m in self.MetaList:
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
                MetaList=self.MetaList,
                Surfaces=self.Surfaces,
                UniverseBox=self.UniverseBox,
                debug=self.settings.debug,
                meta=True,
                tolerances=self.tolerances,
                options=self.options,
                numeric_format=self.numeric_format
            )

            # decompose Enclosure solids
            if self.settings.voidGen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    MetaList=EnclosureList,
                    Surfaces=self.Surfaces,
                    UniverseBox=self.UniverseBox,
                    debug=self.settings.debug,
                    meta=False,
                    tolerances=self.tolerances,
                    options=self.options,
                    numeric_format=self.numeric_format
                )

            if self.options.verbose:
                print("End of decomposition phase")

            # start Building CGS cells phase
            for j, m in enumerate(self.MetaList):
                if m.IsEnclosure:
                    continue
                if self.options.verbose:
                    print("Building cell: ", j + 1)
                cones = Conv.cellDef(
                    meta_obj=m,
                    surfaces=self.Surfaces,
                    universe_box=self.UniverseBox,
                    tolerances=self.tolerances,
                    options=self.options,
                    numeric_format=self.numeric_format
                )
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningSolidList:
                    warnSolids.append(m)
                if not m.Solids:
                    if self.options.verbose:
                        print("none", j, m.__id__)
                        print(m.Definition)

            if self.options.forceNoOverlap:
                Conv.no_overlapping_cell(self.MetaList, self.Surfaces)

        else:
            translate(self.MetaList, self.Surfaces, self.UniverseBox, self.debug, self.options.verbose)
            # decompose Enclosure solids
            if self.voidGen and EnclosureList:
                warningEnclosureList = decompose_solids(
                    MataList=EnclosureList,
                    Surfaces=self.Surfaces,
                    UniverseBox=self.UniverseBox,
                    debug=self.debug,
                    meta=False,
                    tolerances=self.tolerances,
                    options=self.options,
                    numeric_format=self.numeric_format
                )

        if self.options.verbose:
            tempstr2 = str(datetime.now() - tempTime)
            print(f"CSG cells built, time = {tempstr2}")

        #  building enclosure solids
        if self.settings.voidGen and EnclosureList:
            for j, m in enumerate(EnclosureList):
                if self.options.verbose:
                    print("Building Enclosure Cell: ", j + 1)
                cones = Conv.cellDef(
                    meta_obj=m,
                    surfaces=self.Surfaces,
                    universe_box=self.UniverseBox,
                    tolerances=self.tolerances,
                    options=self.options,
                )
                if cones:
                    coneInfo[m.__id__] = cones
                if j in warningEnclosureList:
                    warnEnclosures.append(m)

        tempTime1 = datetime.now()

        # void generation phase
        MetaVoid = []
        if self.settings.voidGen:
            if self.options.verbose:
                print("Build Void")
                print(self.settings.voidExclude)
            if not self.settings.voidExclude:
                MetaReduced = self.MetaList
            else:
                MetaReduced = exclude_cells(self.MetaList, self.settings.voidExclude)

            if self.MetaList:
                init = self.MetaList[-1].__id__ - len(EnclosureList)
            else:
                init = 0
            # TODO perhaps this method should be moved into the CsgToCsg class to avoid passing in so many args
            MetaVoid = Void.void_generation(
                MetaList=MetaReduced,
                EnclosureList=EnclosureList,
                Surfaces=self.Surfaces,
                UniverseBox=self.UniverseBox,
                init=init,
                settings=self.settings,
                tolerances=self.tolerances,
                options=self.options,
                numeric_format=self.numeric_format
            )

        # if self.simplify == 'full' and not Options.forceNoOverlap:
        if self.settings.simplify == "full":
            Surfs = {}
            for lst in self.Surfaces.values():
                for s in lst:
                    Surfs[s.Index] = s

            for c in self.MetaList:
                if c.Definition.level == 0 or c.IsEnclosure:
                    continue
                if self.options.verbose:
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
                    if self.options.verbose:
                        print(f"unexpected constant cell {c.__id__} :{c.Definition.elements}")

        if self.options.verbose:
            tempTime2 = datetime.now()
            print("build Time:", tempTime2 - tempTime1)

            print(datetime.now() - startTime)

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

        if self.options.verbose:
            print_warning_solids(warnSolids, warnEnclosures)

            # add plane definition to cone
            process_cones(self.MetaList, coneInfo, self.Surfaces, self.UniverseBox, self.tolerances.angle)

            print(f"Process finished, time = {datetime.now() - startTime}")

            print("Translation time of solid cells", tempTime1 - tempTime0)
            print("Translation time of void cells", tempTime2 - tempTime1)
   
    def export_csg(
            self,
            title: str = "Geouned conversion",
            geometry_name: str = "converted_with_geouned",
            out_formats: typing.Tuple[str] = ("openmc_xml", "openmc_py", "serpent", "phits", "mcnp"),
            volSDEF: bool = False,
            volCARD: bool = True,
            UCARD=None,
            dummyMat: bool = False,
            cellCommentFile: bool = False,
            cellSummaryFile: bool = True,
    ):
        """Writes out a CSG file in the requested Monte Carlo code format

        Args:
            title (str, optional): Title of the model written at the top of the output file. Defaults to "Geouned conversion".
            geometry_name (str, optional): file stem of the output file(s). Defaults to "converted_with_geouned".
            out_format (typing.Tuple[str], optional): Format for the output
                geometry. Available format are: "mcnp", "openmc_xml",
                "openmc_py", "phits" and "serpent". Several output format can
                be written in the same method call. Defaults to output all codes.
                ("openmc_xml", "openmc_py", "serpent", "phits", "mcnp").
            volSDEF (bool, optional):  Write SDEF definition and tally of solid
                cell for stochastic volume checking. Defaults to False.
            volCARD (bool, optional): Write the CAD calculated volume in the
                cell definition using the VOL card. Defaults to True.
            UCARD (_type_, optional): Write universe card in the cell definition
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

        supported_mc_codes = ("mcnp", "openmc_xml", "openmc_py", "serpent", "phits")
        for out_format in out_formats:
            if out_format not in supported_mc_codes:
                msg = f"outFormat {out_format} not in supported MC codes. The supported codes are {supported_mc_codes}"
                raise ValueError(msg)

        files_written = write_geometry(
            StepFile=self.stepFile,
            volSDEF=volSDEF,
            title=title,
            UCARD=UCARD,
            volCARD=volCARD,
            geometryName=geometry_name,
            outFormat=out_formats,
            dummyMat=dummyMat,
            cellCommentFile=cellCommentFile,
            cellSummaryFile=cellSummaryFile,
            UniverseBox=self.UniverseBox,
            MetaList=self.MetaList,
            Surfaces=self.Surfaces,
            settings=self.settings,
            tolerances=self.tolerances,
            numeric_format=self.numeric_format,
            options=self.options
        )
        print(f"Written CSG geometry files {files_written}")



def decompose_solids(
    MetaList,
    Surfaces,
    UniverseBox,
    debug,
    meta,
    tolerances,
    options,
    numeric_format
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
            universe_box=UniverseBox,
            tolerances=tolerances,
            options=options,
            numeric_format=numeric_format
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
        ext_surfaces = Decom.extract_surfaces(
            solid=comsolid,
            kind="All",
            universe_box=UniverseBox,
            tolerances=tolerances,
            numeric_format=numeric_format,
            options=options,
            MakeObj=True,
        )
        Surfaces.extend(
            surface=ext_surfaces,
            tolerances=tolerances,
            options=options,numeric_format=numeric_format
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


def process_cones(MetaList, coneInfo, Surfaces, UniverseBox, angle):
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
            Conv.add_cone_plane(m.Definition, cones, Surfaces, UniverseBox, angle)
        elif not m.Void:
            Conv.add_cone_plane(m.Definition, coneInfo[m.__id__], Surfaces, UniverseBox, angle)


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