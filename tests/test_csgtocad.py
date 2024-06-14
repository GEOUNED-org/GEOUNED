import geouned

from pathlib import Path

@pytest.mark.parametrize("csg_format", ['mcnp', 'openmc_xml'])
def test_cylbox_convertion(csg_format):
    geo = geouned.CsgToCad()

    # testing/inputSTEP/cylBox.stp has a rough bounding box of
    # -1000.0,  -500.0, -1000.0, 0,0,0.0 in mm
    geo.export_cad(
        input_filename='tests_outputs/testing/inputSTEP/cylBox/cylBox.mcnp',
        csg_format=csg_format,
        bounding_box=[-1000.0,  -500.0, -1000.0, 0,0,0.0 ],
        # cell_range_type='exclude',
        # cell_range=(2,3,4),
        output_filename=f'tests_outputs/csgtocad/cylBox_{csg_format}.stp'
    )
    assert Path('tests_outputs/csgtocad/cylBox.stp').exists()
