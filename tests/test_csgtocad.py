from pathlib import Path
import pytest
import geouned


@pytest.mark.parametrize("csg_format", ["mcnp", "openmc_xml"])
def test_cylbox_convertion(csg_format):

    if csg_format == "openmc_xml":
        suffix = ".xml"
    elif csg_format == "mcnp":
        suffix = ".mcnp"

    geo = geouned.CsgToCad()

    geo.read_csg_file(
        # csg file was made from testing/inputSTEP/cylBox.stp
        input_filename=f"tests/csg_files/cylinder_box{suffix}",
        csg_format=csg_format,
    )

    geo.build_universe()
    geo.export_cad(output_filename=f"tests_outputs/csgtocad/{csg_format}")

    assert Path(f"tests_outputs/csgtocad/{csg_format}.stp").exists()
    assert Path(f"tests_outputs/csgtocad/{csg_format}.FCStd").exists()
