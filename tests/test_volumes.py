"""
The tests check the resulting volume of the CSG cells using OpenMC.
The tests assumes that the test_convert script has previously been run and
that xml files exist in the tests_outputs folder.
"""

import math
import os
import sys
from pathlib import Path

try:
    import freecad  # importing conda package if present
except:
    pass
import geouned
import openmc
import Part
import pytest
from FreeCAD import Import

with open("cross_sections.xml", "w") as file:
    file.write(
        """
        <?xml version='1.0' encoding='UTF-8'?>
        <cross_sections>
        <library materials="H1" path="" type="neutron"/>
        </cross_sections>
    """
    )
openmc.config["cross_sections"] = Path("cross_sections.xml").resolve()

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))
# this file is removed as it fails the test. An issue has been raised to track
step_files.remove(Path("testing/inputSTEP/Misc/rails.stp"))
# this checks if the tests are being run a github action runner or not
if os.getenv("GITHUB_ACTIONS"):
    # reduced samples as github action runners have 2 threads and more samples
    # is not needed for the smaller models tested. This allows the CI to
    # quickly test the smaller models
    samples = 4_000_000
    rel_tol = 0.05
    # removing two geometries that are particularly slow to convert from CI testing
    # these two geometries remain in the test suite for locally testing
    step_files.remove(Path("testing/inputSTEP/large/SCDR.stp"))
    step_files.remove(Path("testing/inputSTEP/large/Triangle.stp"))
else:
    # samples for local run can be larger as threads is likely to be larger
    samples = 40_000_000
    # acceptable tolerance can also be smaller
    rel_tol = 0.04


@pytest.mark.skipif(
    sys.platform in ["win32", "darwin"],
    reason="OpenMC doesn't install on Windows currently and is not well tested on Mac",
)
@pytest.mark.parametrize("input_step_file", step_files)
def test_volumes(input_step_file):

    output_dir = Path("tests_outputs") / input_step_file.with_suffix("")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename_stem = output_dir / input_step_file.stem
    openmc_xml_file = output_filename_stem.with_suffix(".xml")

    Import.insert(str(input_step_file), "converted-cad")
    result = Part.Shape()
    result.read(str(input_step_file))

    materials = openmc.Materials()
    geometry = openmc.Geometry.from_xml(openmc_xml_file, materials)
    cells = geometry.get_all_cells()

    # some of the csg cells are void space around the CAD, so we just check
    # that there are at least as many cells as there could be more cells in the
    # CSG than the CAD
    assert len(result.Solids) <= len(cells.keys())

    all_cell_vol_calc = []
    for solid, (cell_id, cell) in zip(result.Solids, cells.items()):

        bb = solid.BoundBox
        llx, lly, llz, urx, ury, urz = (
            bb.XMin,
            bb.YMin,
            bb.ZMin,
            bb.XMax,
            bb.YMax,
            bb.ZMax,
        )
        llx, lly, llz, urx, ury, urz = (
            llx / 10,
            lly / 10,
            llz / 10,
            urx / 10,
            ury / 10,
            urz / 10,
        )

        cell_vol_calc = openmc.VolumeCalculation(
            domains=[cell],
            samples=samples,
            lower_left=(llx, lly, llz),
            upper_right=(urx, ury, urz),
        )
        # could add a trigger to ends tests early once they have converged
        # cell_vol_calc.set_trigger(0.05, 'rel_err')

        all_cell_vol_calc.append(cell_vol_calc)

    settings = openmc.Settings()
    settings.volume_calculations = all_cell_vol_calc
    settings.run_mode = "volume"
    model = openmc.Model(geometry, materials, settings)

    model.export_to_model_xml(path=openmc_xml_file.with_name("model.xml"))
    openmc.run(cwd=output_dir)

    for solid, (cell_id, cell) in zip(result.Solids, cells.items()):

        vol_calc_file = openmc_xml_file.with_name(f"volume_{cell_id}.h5")
        cell_vol_calc_results = openmc.VolumeCalculation.from_hdf5(vol_calc_file)
        volume_of_csg_cell = cell_vol_calc_results.volumes[cell_id].nominal_value
        # converts from mm3 in cad to cm3 in csg
        volume_of_cad_cell = solid.Volume * 0.001
        assert math.isclose(volume_of_cad_cell, volume_of_csg_cell, rel_tol=rel_tol)


@pytest.mark.parametrize("skip_solids,expected", [([], 12), ([1, 4, 6], 9)])
def test_skip_solids_when_loading(skip_solids, expected):
    """test to check that the correct number of cells are found in the CSG file when skip solids is set"""

    geo = geouned.CadToCsg(settings=geouned.Settings(voidGen=False, compSolids=False))

    # this geometry has 12 solids in it
    geo.load_step_file(filename="testing/inputSTEP/Torus/example.stp", skip_solids=skip_solids)
    geo.start()
    assert geo.filename == "testing/inputSTEP/Torus/example.stp"
    assert geo.skip_solids == skip_solids

    geo.export_csg(
        geometryName="tests_outputs/skip_solids/csg",
        outFormat=["openmc_xml"],
    )

    materials = openmc.Materials()
    geometry = openmc.Geometry.from_xml("tests_outputs/skip_solids/csg.xml", materials)
    cells = geometry.get_all_cells()
    all_cell_ids = cells.keys()
    assert len(all_cell_ids) == expected
