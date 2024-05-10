from . import AdditionalFiles as OutFiles
from .MCNPFormat import McnpInput
from .OpenMCFormat import OpenmcInput
from .PHITSFormat import PhitsInput
from .SerpentFormat import SerpentInput


def write_geometry(
    universe_box,
    meta_list,
    surfaces,
    filename_stem,
    out_formats,
    cell_comment_file,
    cell_summary_file,
    code_setting,
    options,
    tolerances,
    numeric_format,
):

    # write cells comments in file
    if cell_comment_file:
        OutFiles.comments_write(filename_stem, meta_list)
    if cell_summary_file:
        OutFiles.summary_write(filename_stem, meta_list)

    if "mcnp" in out_formats:
        mcnpFilename = filename_stem + ".mcnp"
        outBox = (
            universe_box.XMin,
            universe_box.XMax,
            universe_box.YMin,
            universe_box.YMax,
            universe_box.ZMin,
            universe_box.ZMax,
        )
        if code_setting["voidGen"]:
            outSphere = (surfaces["Sph"][-1].Index, surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        MCNPfile = McnpInput(
            meta_list, surfaces, code_setting, options, tolerances, numeric_format
        )
        MCNPfile.set_sdef((outSphere, outBox))
        MCNPfile.write_input(mcnpFilename)

    if "openMC_XML" in out_formats or "openMC_PY" in out_formats:
        OMCFile = OpenmcInput(meta_list, surfaces, options, tolerances, numeric_format)

    if "openMC_XML" in out_formats:
        omcFilename = filename_stem + ".xml"
        OMCFile.write_xml(omcFilename)

    if "openMC_PY" in out_formats:
        omcFilename = filename_stem + ".py"
        OMCFile.write_py(omcFilename)

    if "serpent" in out_formats:
        serpentFilename = filename_stem + ".serp"
        outBox = (
            universe_box.XMin,
            universe_box.XMax,
            universe_box.YMin,
            universe_box.YMax,
            universe_box.ZMin,
            universe_box.ZMax,
        )
        if code_setting["voidGen"]:
            outSphere = (surfaces["Sph"][-1].Index, surfaces["Sph"][-1].Surf.Radius)
        else:
            outSphere = None

        Serpentfile = SerpentInput(
            meta_list, surfaces, code_setting, options, tolerances, numeric_format
        )
        # Serpentfile.set_sdef((outSphere,outBox))
        Serpentfile.write_input(serpentFilename)

    if "phits" in out_formats:
        phitsFilename = filename_stem + ".inp"
        PHITS_outBox = (
            universe_box.XMin,
            universe_box.XMax,
            universe_box.YMin,
            universe_box.YMax,
            universe_box.ZMin,
            universe_box.ZMax,
        )
        if code_setting["voidGen"]:
            PHITS_outSphere = (
                surfaces["Sph"][-1].Index,
                surfaces["Sph"][-1].Surf.Radius,
            )
        else:
            PHITS_outSphere = None

        PHITSfile = PhitsInput(
            meta_list, surfaces, code_setting, options, tolerances, numeric_format
        )
        # PHITSfile.setSDEF_PHITS((PHITS_outSphere,PHITS_outBox))
        PHITSfile.write_phits(phitsFilename)
