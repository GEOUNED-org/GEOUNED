# Examples

## Conversion of CAD to CSG


The following code snippet will convert a CAD file in the STEP format called example.step into and XML file called geometry.xml that can be opened and run with OpenMC

    import geouned
    geouned.cad_to_csg(
        cad_filename='example.step',
        output_filename='geometry.xml',
        csg_format='openmc'
    )