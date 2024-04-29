from pathlib import Path

import pytest

from geouned import CadToCsg

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))


@pytest.mark.parametrize("input_step_file", step_files)
def test_conversion(input_step_file):
    """Test that step files can be converted to openmc and mcnp files"""

    # sets up an output folder for the results
    output_dir = Path("tests_outputs") / input_step_file.with_suffix("")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename_stem = output_dir / input_step_file.stem

    # creates the config file contents
    template = {
        "title": "Input Test",
        "stepFile": f"{input_step_file.resolve()}",
        "geometryName": f"{output_filename_stem.resolve()}",
        "outFormat": ("mcnp", "openMC_XML"),
        "compSolids": False,
        "volCARD": False,
        "volSDEF": True,
        "voidGen": True,
        "dummyMat": True,
        "minVoidSize": 100,
        "cellSummaryFile": False,
        "cellCommentFile": False,
        "debug": False,
        "simplify": "no",
        "forceCylinder": False,
        "splitTolerance": 0,
        "newSplitPlane": True,
        "nPlaneReverse": 0,
    }

    # deletes the output openmc and mcnp output files if it already exists
    output_filename_stem.with_suffix(".mcnp").unlink(missing_ok=True)
    output_filename_stem.with_suffix(".xml").unlink(missing_ok=True)

    GEO = CadToCsg("Input Test")

    # set parameters values stored in template dictionary
    for key, value in template.items():
        GEO.set(key, value)

    GEO.Start()

    assert output_filename_stem.with_suffix(".mcnp").exists()
    assert output_filename_stem.with_suffix(".xml").exists()
