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

    geo.export_cad(
        # csg file was made from testing/inputSTEP/cylBox.stp
        input_filename=f"tests/csg_files/cylinder_box{suffix}",
        csg_format=csg_format,
        bounding_box=[-1000.0, -500.0, -1000.0, 0, 0, 0.0],
        # TODO add tests for these args that counts volumes in cad file
        # cell_range_type='exclude',
        # cell_range=(2,3,4),
        output_filename=f"tests_outputs/csgtocad/{csg_format}",
    )

    assert Path(f"tests_outputs/csgtocad/{csg_format}.step").exists()
    assert Path(f"tests_outputs/csgtocad/{csg_format}.FCStd").exists()
