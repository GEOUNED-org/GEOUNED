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
    
    my_settings = geouned.Settings(
        compSolids=False,
        minVoidSize=100,
        debug=False,
        simplify="no",
    )

    geo = geouned.CadToCsg(
        stepFile=f"{input_step_file.resolve()}",
        options=my_options,
        settings=my_settings,
    )

    geo.start()
    geo.export_csg(
        cellSummaryFile=False,
        cellCommentFile=False,
        dummyMat=True,
        title="Input Test",
        geometry_name=f"{output_filename_stem.resolve()}",
        out_formats=("mcnp", "openmc_xml", "openmc_py", "serpent", "phits"),
        volCARD=False,
        volSDEF=True,
    )

    for suffix in suffixes:
        assert output_filename_stem.with_suffix(suffix).exists()
