from . import AdditionalFiles as OutFiles
from .MCNPFormat import McnpInput
from .OpenMCFormat import OpenmcInput
from .PHITSFormat import PhitsInput
from .SerpentFormat import SerpentInput


def write_geometry(UniverseBox, MetaList, Surfaces, code_setting):

    baseName = code_setting["geometryName"]

    # Currently there are two was of setting outFormat (via a .set method and
    # a class attribute. Once we have a single method then move this validating
    # input code to the attribute @setter
    supported_mc_codes = ("mcnp", "openMC_XML", "openMC_PY", "serpent", "phits")
    for out_format in code_setting["outFormat"]:
        if out_format not in supported_mc_codes:
            msg = f"outFormat {out_format} not in supported MC codes ({supported_mc_codes})"
            raise ValueError(msg)

    # write cells comments in file
    if code_setting["cellCommentFile"]:
        OutFiles.comments_write(baseName, MetaList)
    if code_setting["cellSummaryFile"]:
        OutFiles.summary_write(baseName, MetaList)

    if "mcnp" in code_setting["outFormat"]:
        mcnpFilename = baseName + ".mcnp"
        outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if code_setting["voidGen"]:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        MCNPfile = McnpInput(MetaList, Surfaces, code_setting)
        MCNPfile.set_sdef((outSphere, outBox))
        MCNPfile.write_input(mcnpFilename)

    if (
        "openMC_XML" in code_setting["outFormat"]
        or "openMC_PY" in code_setting["outFormat"]
    ):
        OMCFile = OpenmcInput(MetaList, Surfaces, code_setting)

    if "openMC_XML" in code_setting["outFormat"]:
        omcFilename = baseName + ".xml"
        OMCFile.write_xml(omcFilename)

    if "openMC_PY" in code_setting["outFormat"]:
        omcFilename = baseName + ".py"
        OMCFile.write_py(omcFilename)

    if "serpent" in code_setting["outFormat"]:
        serpentFilename = baseName + ".serp"
        outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if code_setting["voidGen"]:
            outSphere = (Surfaces["Sph"][-1].Index, Surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        Serpentfile = SerpentInput(MetaList, Surfaces, code_setting)
        # Serpentfile.set_sdef((outSphere,outBox))
        Serpentfile.write_input(serpentFilename)

    if "phits" in code_setting["outFormat"]:
        phitsFilename = baseName + ".inp"
        PHITS_outBox = (
            UniverseBox.XMin,
            UniverseBox.XMax,
            UniverseBox.YMin,
            UniverseBox.YMax,
            UniverseBox.ZMin,
            UniverseBox.ZMax,
        )
        if code_setting["voidGen"]:
            PHITS_outSphere = (
                Surfaces["Sph"][-1].Index,
                Surfaces["Sph"][-1].Surf.Radius,
            )
        else:
            PHITS_outSphere = None

        PHITSfile = PhitsInput(MetaList, Surfaces, code_setting)
        # PHITSfile.setSDEF_PHITS((PHITS_outSphere,PHITS_outBox))
        PHITSfile.write_phits(phitsFilename)
