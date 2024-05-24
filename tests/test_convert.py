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
suffixes = (".mcnp", ".xml", ".inp", ".py", ".serp")


@pytest.mark.parametrize("input_step_file", step_files)
def test_conversion(input_step_file):
    """Test that step files can be converted to openmc and mcnp files"""

    # sets up an output folder for the results
    output_dir = Path("tests_outputs") / input_step_file.with_suffix("")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename_stem = output_dir / input_step_file.stem

    # deletes the output MC files if they already exists
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

    my_settings = geouned.Settings(
        matFile="",
        voidGen=True,
        debug=False,
        compSolids=True,
        simplify="no",
        cellRange=[],
        exportSolids="",
        minVoidSize=200.0,  # units mm
        maxSurf=50,
        maxBracket=30,
        voidMat=[],
        voidExclude=[],
        startCell=1,
        startSurf=1,
        sort_enclosure=False,
    )

    geo = geouned.CadToCsg(
        stepFile=f"{input_step_file.resolve()}",
        options=my_options,
        settings=my_settings,
        tolerances=my_tolerances,
        numeric_format=my_numeric_format,
    )

    geo.start()

    geo.export_csg(
        title="Converted with GEOUNED",
        geometryName=f"{output_filename_stem.resolve()}",
        outFormat=(
            "openmc_xml",
            "openmc_py",
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

    for suffix in suffixes:
        assert output_filename_stem.with_suffix(suffix).exists()


@pytest.mark.parametrize(
    "input_json_file",
    ["tests/config_complete_defaults.json", "tests/config_minimal.json"],
)
def test_cad_to_csg_from_json_with_defaults(input_json_file):

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        Path("csg").with_suffix(suffix).unlink(missing_ok=True)

    my_cad_to_csg = geouned.CadToCsg.from_json(input_json_file)
    assert isinstance(my_cad_to_csg, geouned.CadToCsg)

    assert my_cad_to_csg.stepFile == "testing/inputSTEP/BC.stp"
    assert my_cad_to_csg.options.forceCylinder == False
    assert my_cad_to_csg.tolerances.relativeTol == False
    assert my_cad_to_csg.numeric_format.P_abc == "14.7e"
    assert my_cad_to_csg.settings.matFile == ""

    for suffix in suffixes:
        assert Path("csg").with_suffix(suffix).exists()

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        Path("csg").with_suffix(suffix).unlink(missing_ok=True)

    my_cad_to_csg.start()
    my_cad_to_csg.export_csg()


def test_cad_to_csg_from_json_with_non_defaults():

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        Path("csg").with_suffix(suffix).unlink(missing_ok=True)

    my_cad_to_csg = geouned.CadToCsg.from_json("tests/config_non_defaults.json")
    assert isinstance(my_cad_to_csg, geouned.CadToCsg)

    assert my_cad_to_csg.stepFile == "testing/inputSTEP/BC.stp"
    assert my_cad_to_csg.options.forceCylinder == True
    assert my_cad_to_csg.tolerances.relativePrecision == 2e-6
    assert my_cad_to_csg.numeric_format.P_abc == "15.7e"
    assert my_cad_to_csg.settings.matFile == "non default"

    for suffix in suffixes:
        assert Path("csg").with_suffix(suffix).exists()

    # deletes the output MC files if they already exists
    for suffix in suffixes:
        Path("csg").with_suffix(suffix).unlink(missing_ok=True)

    my_cad_to_csg.start()
    my_cad_to_csg.export_csg()


def test_writing_to_new_folders():
    """Checks that a folder is created prior to writing output files"""

    geo = geouned.CadToCsg(stepFile="testing/inputSTEP/BC.stp")
    geo.start()

    for outformat in ["mcnp", "phits", "serpent", "openmc_xml", "openmc_py"]:
        geo.export_csg(
            geometryName=f"tests_outputs/new_folder_for_testing_{outformat}/csg",
            cellCommentFile=False,
            cellSummaryFile=False,
            outFormat=[outformat],
        )
        geo.export_csg(
            geometryName=f"tests_outputs/new_folder_for_testing_{outformat}_cell_comment/csg",
            cellCommentFile=True,
            cellSummaryFile=False,
            outFormat=[outformat],
        )
        geo.export_csg(
            geometryName=f"tests_outputs/new_folder_for_testing_{outformat}_cell_summary/csg",
            cellCommentFile=False,
            cellSummaryFile=True,
            outFormat=[outformat],
        )


def test_with_relative_tol_true():

    # test to protect against incorrect attribute usage in FreeCAD
    # more details https://github.com/GEOUNED-org/GEOUNED/issues/154

    geo = geouned.CadToCsg(
        stepFile=f"{step_files[1].resolve()}",
        tolerances=geouned.Tolerances(relativeTol=False),
    )
    geo.start()
    geo = geouned.CadToCsg(
        stepFile=f"{step_files[1].resolve()}",
        tolerances=geouned.Tolerances(relativeTol=True),
    )
    geo.start()


@pytest.mark.parametrize("input_step_file", step_files)
@pytest.mark.parametrize("suffix", suffixes)
def test_new_mc_files_match_original(suffix, input_step_file):
    """
    Regression test that the checks the text files produced for each MC code match the text files produced previously.
    If this test fails it might be due to an improved MC code output instead of a mistake in the PR.
    You might want to update the MC text file in the regression test folder with the 'tests/update_regression_test_files.py' script.
    """

    # sets up an output folder for the results
    regression_test_file = (
        Path("tests/regression_test_files")
        / input_step_file.parts[-2]
        / Path(input_step_file.stem)
        / Path(input_step_file.name).with_suffix(suffix)
    )
    output_filename = Path("tests_outputs") / input_step_file.with_suffix("") / Path(input_step_file.stem).with_suffix(suffix)

    with open(output_filename, "r") as f:
        file_new = f.readlines()
    with open(regression_test_file, "r") as f:
        file_original = f.readlines()
    for line_new, line_original in zip(file_new, file_original):
        # this lines in the output files are not expected to match
        if (
            " Creation Date" not in line_new
            and " Creation Date" not in line_original
            and " Version : " not in line_new
            and " Version : " not in line_original
            and " Original Step file :  " not in line_new
            and " Original Step file :  " not in line_original
        ):
            assert line_new == line_original
    assert len(file_new) == len(file_original)
