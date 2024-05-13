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
        newSplitPlane=True,
        delLastNumber=False,
        enlargeBox=2,
        nPlaneReverse=0,
        splitTolerance=0,
        scaleUp=True,
        quadricPY=False,
        Facets=False,
        prnt3PPlane=False,
        forceNoOverlap=False,
    )

    my_tolerances = geouned.Tolerances(
        relativeTol=False,
        relativePrecision=0.000001,
        value=0.000001,
        distance=0.0001,
        angle=0.0001,
        pln_distance=0.0001,
        pln_angle=0.0001,
        cyl_distance=0.0001,
        cyl_angle=0.0001,
        sph_distance=0.0001,
        kne_distance=0.0001,
        kne_angle=0.0001,
        tor_distance=0.0001,
        tor_angle=0.0001,
        min_area=0.01,
    )
    my_numeric_format = geouned.NumericFormat(
        P_abc="14.7e",
        P_d="14.7e",
        P_xyz="14.7e",
        S_r="14.7e",
        S_xyz="14.7e",
        C_r="12f",
        C_xyz="12f",
        K_xyz="13.6e",
        K_tan2="12f",
        T_r="14.7e",
        T_xyz="14.7e",
        GQ_1to6="18.15f",
        GQ_7to9="18.15f",
        GQ_10="18.15f",
    )

    geo = geouned.CadToCsg(
        stepFile=f"{input_step_file.resolve()}",
        matFile="",
        voidGen=True,
        debug=False,
        compSolids=False,  # changed from the default
        simplify="no",
        cellRange=[],
        exportSolids="",
        minVoidSize=100,  # changed from the default
        maxSurf=50,
        maxBracket=30,
        voidMat=[],
        voidExclude=[],
        startCell=1,
        startSurf=1,
        sort_enclosure=False,
        options=my_options,
        tolerances=my_tolerances,
        numeric_format=my_numeric_format,
        title="Converted with GEOUNED",
        geometryName=f"{output_filename_stem.resolve()}",
        outFormat=(
            "openMC_XML",
            "openMC_PY",
            "serpent",
            "phits",
            "mcnp",
        ),
        volSDEF=True,  # changed from the default
        volCARD=False,  # changed from the default
        UCARD=None,
        dummyMat=True,  # changed from the default
        cellCommentFile=False,
        cellSummaryFile=False,  # changed from the default
    )

    geo.start()

    for suffix in suffixes:
        assert output_filename_stem.with_suffix(suffix).exists()
