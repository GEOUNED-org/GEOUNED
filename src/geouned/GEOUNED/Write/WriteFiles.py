from . import AdditionalFiles as OutFiles
from .MCNPFormat import McnpInput
from .OpenMCFormat import OpenmcInput
from .PHITSFormat import PhitsInput
from .SerpentFormat import SerpentInput


def write_geometry(
    step_file,
    volSDEF,
    title,
    UCARD,
    volCARD,
    geometryName,
    outFormat,
    dummyMat,
    cellCommentFile,
    cellSummaryFile,
    UniverseBox,
    MetaList,
    Surfaces,
    settings,
    tolerances,
    numeric_format,
    options,
):

    files_written = []

    # write cells comments in file
    if cellCommentFile:
        OutFiles.comments_write(geometryName, MetaList)
    if cellSummaryFile:
        OutFiles.summary_write(geometryName, MetaList)

    if "mcnp" in outFormat:
        mcnpFilename = geometryName + ".mcnp"
        outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if settings.voidGen:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        MCNPfile = McnpInput(
            MetaList,
            Surfaces,
            title,
            volSDEF,
            volCARD,
            UCARD,
            dummyMat,
            options.prnt3PPlane,
            step_file,
            tolerances,
            numeric_format,
            options=options,
        )
        MCNPfile.set_sdef((outSphere, outBox))
        MCNPfile.write_input(mcnpFilename, options.verbose)
        files_written.append(mcnpFilename)

    if "openmc_xml" in outFormat or "openmc_py" in outFormat:
        OMCFile = OpenmcInput(
            Meta=MetaList,
            Surfaces=Surfaces,
            tolerances=tolerances,
            numeric_format=numeric_format,
            options=options,
        )

    if "openmc_xml" in outFormat:
        omcFilename = geometryName + ".xml"
        OMCFile.write_xml(omcFilename, options.verbose)
        files_written.append(omcFilename)

    if "openmc_py" in outFormat:
        omcFilename = geometryName + ".py"
        OMCFile.write_py(omcFilename, options.verbose)
        files_written.append(omcFilename)

    if "serpent" in outFormat:
        serpentFilename = geometryName + ".serp"
        outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if settings.voidGen:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        Serpentfile = SerpentInput(
            MetaList,
            Surfaces,
            title,
            volSDEF,
            volCARD,
            UCARD,
            dummyMat,
            options.prnt3PPlane,
            step_file,
            tolerances,
            numeric_format,
            options,
        )
        # Serpentfile.set_sdef((outSphere,outBox))
        Serpentfile.write_input(serpentFilename, options.verbose)
        files_written.append(serpentFilename)

    if "phits" in outFormat:
        phitsFilename = geometryName + ".inp"
        PHITS_outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if settings.voidGen:
            PHITS_outSphere = (
                Surfaces["Sph"][-1].Index,
                Surfaces["Sph"][-1].Surf.Radius,
            )
        else:
            PHITS_outSphere = None

        PHITSfile = PhitsInput(
            MetaList,
            Surfaces,
            title,
            volSDEF,
            volCARD,
            UCARD,
            dummyMat,
            settings.matFile,
            settings.voidMat,
            settings.startCell,
            options.prnt3PPlane,
            step_file,
            tolerances,
            numeric_format,
            options,
        )
        # PHITSfile.setSDEF_PHITS((PHITS_outSphere,PHITS_outBox))
        PHITSfile.write_phits(phitsFilename, options.verbose)
        files_written.append(phitsFilename)
    return files_written
