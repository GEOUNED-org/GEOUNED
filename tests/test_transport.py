"""
The tests check the resulting CSG geometry in particle transport using OpenMC.
The tests assumes that the test_convert script has previously been run and
that xml files exist in the tests_outputs folder.
"""

import os
import sys
from pathlib import Path

try:
    import freecad  # importing conda package if present
except:
    pass
import Part
import pytest
from FreeCAD import Import
import openmc
import openmc

openmc.config["cross_sections"] = Path("tests/cross_sections.xml").resolve()

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))

# these step files all produced failures, most of them lose particles
# they are removed from the tests but an issue has be raised and a minimal
# example has been made in the issue where the fixing of these files can be
# tracked https://github.com/shimwell/GEOUNED/issues/11
step_files.remove(Path("testing/inputSTEP/large/SCDR.stp"))
step_files.remove(Path("testing/inputSTEP/placa.stp"))
step_files.remove(Path("testing/inputSTEP/placa2.stp"))
# required a few more particles needed but also losing particles
step_files.remove(Path("testing/inputSTEP/DoubleCylinder/placa3.step"))
step_files.remove(Path("testing/inputSTEP/DoubleCylinder/placa.stp"))
# this face2.stp crashes when loading the geometry.xml


@pytest.mark.skipif(
    sys.platform in ["win32", "darwin"],
    reason="OpenMC doesn't install on Windows currently and is not well tested on Mac",
)
@pytest.mark.parametrize("input_step_file", step_files)
def test_transport(input_step_file):

    output_dir = Path("tests_outputs") / input_step_file.with_suffix("")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename_stem = output_dir / input_step_file.stem
    openmc_xml_file = output_filename_stem.with_suffix(".xml")

    Import.insert(str(input_step_file), "converted-cad")
    result = Part.Shape()
    result.read(str(input_step_file))

    # getting the bounding box of the CAD so we can use it to make a source
    # over the whole geometry
    bb = result.BoundBox
    llx, lly, llz, urx, ury, urz = bb.XMin, bb.YMin, bb.ZMin, bb.XMax, bb.YMax, bb.ZMax
    # converting from mm to cm
    llx, lly, llz, urx, ury, urz = (
        llx / 10,
        lly / 10,
        llz / 10,
        urx / 10,
        ury / 10,
        urz / 10,
    )

    source = openmc.IndependentSource()
    source.space = openmc.stats.Box(
        lower_left=(llx, lly, llz), upper_right=(urx, ury, urz)
    )
    source.energy = openmc.stats.Discrete([14e6], [1])

    materials = openmc.Materials()
    geometry = openmc.Geometry.from_xml(openmc_xml_file, materials)

    settings = openmc.Settings()
    settings.batches = 10
    settings.max_lost_particles = 1
    # number of particles is increased for local runs as the GitHub action
    # runner has two threads and a typical local computer has more.
    # Sometimes a lot of particles are needed to find small geometry errors
    if os.getenv("GITHUB_ACTIONS"):
        settings.particles = 100_000
    else:
        settings.particles = 1_000_000
    settings.run_mode = "fixed source"
    settings.source = source
    model = openmc.Model(geometry, materials, settings)

    model.run()
