Python Package Usage
====================

The main class is ``CsgToCad()`` which converts Constructive Solid Geometry (CSG) to CAD. Class has changed singificantly since version 1.5.x. 
The most relevant modification are: no need to specify external boundbox dimensions anymore, new methods to load, process and write geometry have been added.
The new available methods are ``CsgToCad().read_csg_file()``, ``CsgToCad().build_universe()``, ``CsgToCad().build_container()`` and ``CsgToCad().export_cad()``.
Also methods for cell and material selection were added ``CsgToCad().cell_filter()``, ``CsgToCad().material_filter()``.
The arguments passed into these methods are documented in the `Python API reference section <../python_api.html>`_ of the documentation.

If you have install GEOUNED and FreeCAD into your system Python then you can simply run a .py script with Python.
The most minimal use case below shows GEOUNED being imported and the CsgToCad being used to convert a CSG geometry into a STEP CAD file. 
The example makes use of default attributes.

.. code-block:: python

    import geouned

    geo = geouned.CsgToCad()

    geo.read_csg_file(
        input_filename='cylinder_box.xml',
        csg_format='openmc_xml',
        )
    geo.build_universe()
    geo.export_cad()

Users can change :meth:`geouned.BoxSettings` to modify internal parameters in order to optimize the conversion process.
The following example shows a usage of specidied attributes.
Remember that the arguments are described in the `Python API reference section <../python_api.html>`_ of the documentation.

.. code-block:: python

   import geouned

    my_settings = geouned.BoxSettings(
        universe_radius=1.0e8, 
        insolid_tolerance=0.1, 
    )

    geo = geouned.CsgToCad(settings=my_settings)
    
    geo.read_csg_file(
        input_filename='cylinder_box.xml',
        csg_format='openmc_xml',
        )

    geo.material_filter(
        type='exclude', 
        materials=(0,),
        )
    
    geo.build_universe(
        U=10, 
        depth=-1,
        )

    geo.cell_filter(
        type='include', 
        cells=(101,102,103),
        )
    
    geo.build_container(
        cell_label=100,
        depth=1,
        )

    geo.export_cad(output_filename='cylinder_box.stp')


``CsgToCad().build_universe()`` and ``CsgToCad().build_container()`` are two methods to build the CAD model. 
``CsgToCad().build_universe()`` will build the standalone universe with its own coordinates. Default to root universe. 
If a universe number is passed as argument this universe will be build. The depth argument indicate how deep the nested universe should be considered.
Default to -1 (all nested universes). If a universe is not build because it is in a lower level, its container cell will be included instead.   

``CsgToCad().build_container()`` will build the universe contained the container cell "cell_label". The argurment passed to this method is the cell label in which
the universe is located. This method will build the universe and locate it inside the container with the corresponding transformation (universe region will be 
delimited by the container cell boundaries). In this method the argument "depth" can be used to limit the number of nested universes.