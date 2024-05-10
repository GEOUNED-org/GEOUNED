import os
from pathlib import Path

import pytest

import geouned

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))
# removing two geometries that are particularly slow to convert from CI testing
# these two geometries remain in the test suite for locally testing
if os.getenv("GITHUB_ACTIONS"):
    step_files.remove(Path("testing/inputSTEP/large/SCDR.stp"))
    step_files.remove(Path("testing/inputSTEP/large/Triangle.stp"))


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
        "outFormat": ("mcnp", "openMC_XML", "openMC_PY", "serpent", "phits"),
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
    }

    # deletes the output MC files if they already exists
    suffixes = (".mcnp", ".xml", ".inp", ".py", ".serp")
    for suffix in suffixes:
        output_filename_stem.with_suffix(suffix).unlink(missing_ok=True)

    my_options = geouned.Options(
        forceCylinder=False,
        splitTolerance=0,
        newSplitPlane=True,
        nPlaneReverse=0,
    )

    # default values used, just checking it can be passed in
    my_tolerances = geouned.Tolerances(
        min_area=0.01
    )

    # default values used, just checking it can be passed in
    my_numeric_format = geouned.NumericFormat(
        C_r="12f"
    )

    geo = geouned.CadToCsg(
        title="Input Test", options=my_options, tolerances=my_tolerances, numeric_format=my_numeric_format
    )

    # set parameters values stored in template dictionary
    for key, value in template.items():
        geo.set(key, value)

    geo.start()

    for suffix in suffixes:
        assert output_filename_stem.with_suffix(suffix).exists()
