import json
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
        voidGen=True,
    )

    geo = geouned.CadToCsg(
        step_file=f"{input_step_file.resolve()}",
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


def test_attribute_setting_from_config():

    config_dict = {"step_file": str(step_files[0])}
    with open("config.json", "w") as fp:
        json.dump(config_dict, fp)

    default_geo = geouned.CadToCsg.from_config("config.json")

    assert default_geo.step_file == str(step_files[0])
    assert default_geo.options.forceCylinder == False
    assert default_geo.options.newSplitPlane == True
    assert default_geo.options.delLastNumber == False
    assert default_geo.options.verbose == False
    assert default_geo.options.enlargeBox == 2
    assert default_geo.options.nPlaneReverse == 0
    assert default_geo.options.splitTolerance == 0
    assert default_geo.options.scale_up == True
    assert default_geo.options.quadricPY == False
    assert default_geo.options.Facets == False
    assert default_geo.options.prnt3PPlane == False
    assert default_geo.options.forceNoOverlap == False

    assert default_geo.tolerances.relativeTol == False
    assert default_geo.tolerances.relativePrecision == 1.0e-6
    assert default_geo.tolerances.value == 1.0e-6
    assert default_geo.tolerances.distance == 1.0e-4
    assert default_geo.tolerances.angle == 1.0e-4
    assert default_geo.tolerances.pln_distance == 1.0e-4
    assert default_geo.tolerances.pln_angle == 1.0e-4
    assert default_geo.tolerances.cyl_distance == 1.0e-4
    assert default_geo.tolerances.cyl_angle == 1.0e-4
    assert default_geo.tolerances.sph_distance == 1.0e-4
    assert default_geo.tolerances.kne_distance == 1.0e-4
    assert default_geo.tolerances.kne_angle == 1.0e-4
    assert default_geo.tolerances.tor_distance == 1.0e-4
    assert default_geo.tolerances.tor_angle == 1.0e-4
    assert default_geo.tolerances.min_area == 1.0e-2

    assert default_geo.numeric_format.P_abc == "14.7e"
    assert default_geo.numeric_format.P_d == "14.7e"
    assert default_geo.numeric_format.P_xyz == "14.7e"
    assert default_geo.numeric_format.S_r == "14.7e"
    assert default_geo.numeric_format.S_xyz == "14.7e"
    assert default_geo.numeric_format.C_r == "12f"
    assert default_geo.numeric_format.C_xyz == "12f"
    assert default_geo.numeric_format.K_xyz == "13.6e"
    assert default_geo.numeric_format.K_tan2 == "12f"
    assert default_geo.numeric_format.T_r == "14.7e"
    assert default_geo.numeric_format.T_xyz == "14.7e"
    assert default_geo.numeric_format.GQ_1to6 == "18.15f"
    assert default_geo.numeric_format.GQ_7to9 == "18.15f"
    assert default_geo.numeric_format.GQ_10 == "18.15f"

    assert default_geo.settings.matFile == ""
    assert default_geo.settings.voidGen == True
    assert default_geo.settings.debug == False
    assert default_geo.settings.compSolids == True
    assert default_geo.settings.simplify == "No"
    assert default_geo.settings.cellRange == []
    assert default_geo.settings.exportSolids == ""
    assert default_geo.settings.minVoidSize == 200.0
    assert default_geo.settings.max_surf == 50
    assert default_geo.settings.max_bracket == 30
    assert default_geo.settings.voidMat == []
    assert default_geo.settings.voidExclude == []
    assert default_geo.settings.startCell == 1
    assert default_geo.settings.startSurf == 1
    assert default_geo.settings.sort_enclosure == False

    config_dict = {
        "step_file": str(step_files[0]),
        "Options": {"forceCylinder": True},
        "Tolerances": {"relativePrecision": 1.23e-6},
        "NumericFormat": {"P_abc": "15.6e"},
        "Settings": {"debug": True},
    }
    with open("config.json", "w") as fp:
        json.dump(config_dict, fp, indent=2)

    geo = geouned.CadToCsg.from_config("config.json")
    assert geo.options.forceCylinder == True
    assert geo.tolerances.relativePrecision == 1.23e-6
    assert geo.numeric_format.P_abc == "15.6e"
    assert geo.settings.debug == True
