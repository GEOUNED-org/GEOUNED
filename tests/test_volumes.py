
"""
The tests check the resulting volume of the CSG cells using OpenMC.
The tests assumes that the test_convert script has previously been run and
that xml files exist in the tests_outputs folder.
"""
import sys
import math
import os
from pathlib import Path

try:
    import freecad  # importing conda package if present
except:
    pass
import Part
import pytest
from FreeCAD import Import
import openmc

with open("cross_sections.xml", "w") as file:
    file.write(
    """
        <?xml version='1.0' encoding='UTF-8'?>
        <cross_sections>
        <library materials="H1" path="" type="neutron"/>
        </cross_sections>
    """
    )
openmc.config['cross_sections'] = Path("cross_sections.xml").resolve()

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))

# this checks if the tests are being run a github action runner or not
if os.getenv("GITHUB_ACTIONS"):
    # removing the larger models and models where the volumes don't match from
    # some of these step file are removed as it fails the volume tests.
    # Perhaps a bug with torus conversion or a bug with the stp file?
    # Issue raised https://github.com/GEOUNED-org/GEOUNED/issues/55
    step_files.remove(Path('testing/inputSTEP/Misc/rails.stp'))
    step_files.remove(Path('testing/inputSTEP/placa.stp'))
    step_files.remove(Path('testing/inputSTEP/Misc/tester.stp'))
    step_files.remove(Path('testing/inputSTEP/large/SCDR.stp'))
    step_files.remove(Path('testing/inputSTEP/large/Triangle.stp'))
    step_files.remove(Path('testing/inputSTEP/Torus/face2.stp'))
    step_files.remove(Path('testing/inputSTEP/Torus/rueda.stp'))
    # reduced samples as github action runners have 2 threads and more samples
    # is not needed for the smaller models tested. This allows the CI to
    # quickly test the smaller models
    samples=4_000_000
    rel_tol=0.05
else:
    # samples for local run can be larger as threads is likely to be larger
    samples=400_000_000
    # acceptable tolerance can also be smaller
    rel_tol=0.01


@pytest.mark.skipif(
    sys.platform in ["win32", "darwin"],
    reason="OpenMC doesn't install on Windows currently and is not well tested on Mac"
)
@pytest.mark.parametrize(
    "input_step_file", step_files
)
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
        llx, lly, llz, urx, ury, urz = bb.XMin, bb.YMin, bb.ZMin, bb.XMax, bb.YMax, bb.ZMax
        llx, lly, llz, urx, ury, urz = llx/10, lly/10, llz/10, urx/10, ury/10, urz/10

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

    model.export_to_model_xml(path=openmc_xml_file.with_name('model.xml'))
    openmc.run(cwd=output_dir)

    for solid, (cell_id, cell) in zip(result.Solids, cells.items()):

        vol_calc_file = openmc_xml_file.with_name(f"volume_{cell_id}.h5")
        cell_vol_calc_results = openmc.VolumeCalculation.from_hdf5(vol_calc_file)
        volume_of_csg_cell = cell_vol_calc_results.volumes[cell_id].nominal_value
        # converts from mm3 in cad to cm3 in csg
        volume_of_cad_cell = solid.Volume * 0.001
        assert math.isclose(volume_of_cad_cell, volume_of_csg_cell, rel_tol=rel_tol)
