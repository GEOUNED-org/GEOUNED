Python Package Usage
====================

The main class is ``CsgToCad()`` which converts Constructive Solid Geometry (CSG) to CAD.
There are a few arguments that can be passed into the ``CsgToCad().export_cad()`` method which are documented in the `Python API reference section <../python_api.html>`_ of the documentation.


If you have install GEOUNED and FreeCAD into your system Python then you can simply run a .py script with Python.
The most minimal use case below shows GEOUNED being imported and the CsgToCad being used to convert a CSG geometry into a STEP CAD file. 
The example makes use of default  attributes.

.. code-block:: python

    import geouned

    geo = geouned.CsgToCad()

    geo.export_cad(
        csg_format='openmc_xml',
        input_filename='cylinder_box.xml',
    )


Users can change the default arguments to suit the conversion desired.
The following example shows a usage with every attributes specified.
Remember that the arguments are described in the `Python API reference section <../python_api.html>`_ of the documentation.

.. code-block:: python

    import geouned

    geo = geouned.CsgToCad()

    geo.export_cad(
        input_filename='cylinder_box.xml',
        csg_format='openmc_xml',
        bounding_box=[-1000.0,  -500.0, -1000.0, 0,0,0.0 ],
        cell_range_type='exclude',
        cell_range=(2,3,4),
        output_filename='openmc_xml',
    )
