from . import AdditionalFiles as OutFiles
from .MCNPFormat import McnpInput
from .OpenMCFormat import OpenmcInput
from .PHITSFormat import PhitsInput
from .SerpentFormat import SerpentInput


def write_geometry(
    UniverseBox,
    MetaList,
    Surfaces,
    title,
    volSDEF,
    volCARD,
    UCARD,
    dummyMat,
    geometryName,
    outFormat,
    cellCommentFile,
    cellSummaryFile,
    voidGen,
    matFile,
    voidMat,
    startCell,
):

    # Currently there are two was of setting outFormat (via a .set method and
    # a class attribute. Once we have a single method then move this validating
    # input code to the attribute @setter
    supported_mc_codes = ("mcnp", "openMC_XML", "openMC_PY", "serpent", "phits")
    for out_format in outFormat:
        if out_format not in supported_mc_codes:
            msg = f"outFormat {out_format} not in supported MC codes ({supported_mc_codes})"
            raise ValueError(msg)

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
        if voidGen:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        MCNPfile = McnpInput(
            MetaList, Surfaces, title, volSDEF, volCARD, UCARD, dummyMat
        )
        MCNPfile.set_sdef((outSphere, outBox))
        MCNPfile.write_input(mcnpFilename)

    if "openMC_XML" in outFormat or "openMC_PY" in outFormat:
        OMCFile = OpenmcInput(MetaList, Surfaces)

    if "openMC_XML" in outFormat:
        omcFilename = geometryName + ".xml"
        OMCFile.write_xml(omcFilename)

    if "openMC_PY" in outFormat:
        omcFilename = geometryName + ".py"
        OMCFile.write_py(omcFilename)

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
        if voidGen:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        Serpentfile = SerpentInput(
            MetaList, Surfaces, title, volSDEF, volCARD, UCARD, dummyMat
        )
        # Serpentfile.set_sdef((outSphere,outBox))
        Serpentfile.write_input(serpentFilename)

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
        if voidGen:
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
            matFile,
            voidMat,
            startCell,
        )
        # PHITSfile.setSDEF_PHITS((PHITS_outSphere,PHITS_outBox))
        PHITSfile.write_phits(phitsFilename)
