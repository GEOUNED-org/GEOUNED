Python Package Usage
====================

The main class is ``CadToCsg()`` which converts CAD geometry to Constructive Solid Geometry (CSG).
There are many arguments that can be passed into the ``CadToCsg()`` class which are documented in the `Python API reference section <../python_api.html>`_ of the documentation.


If you have install GEOUNED and FreeCAD into your system Python then you can simply run a .py script with Python.
The most minimal use case below shows GEOUNED being imported and the CadToCsg being used to convert a STEP CAD file called 'cuboid.stp' into a vanity of CSG format. 
The example makes use of default  attributes.

.. code-block:: python

    import geouned
    geo = geouned.CadToCsg()
    geo.load_step_file(filename='cuboid.stp')
    geo.start()
    geo.export_csg()

Users can change :meth:`geouned.Options`, :meth:`geouned.Settings`, :meth:`geouned.Tolerances` and :meth:`geouned.NumericFormat` to suit the conversion desired.
The following example shows a usage with every attributes specified.

.. code-block:: python

    import geouned

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

    my_settings = geouned.Settings(
        matFile="",
        voidGen=True,
        debug=False,
        compSolids=False,
        simplify="no",
        exportSolids="",
        minVoidSize=200.0,
        maxSurf=50,
        maxBracket=30,
        voidMat=[],
        voidExclude=[],
        startCell=1,
        startSurf=1,
        sort_enclosure=False,
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
        options=my_options,
        settings=my_settings,
        tolerances=my_tolerances,
        numeric_format=my_numeric_format,
    )

    geo.load_step_file(
        filename="cuboid.stp",
        skip_solids=[],
    )

    geo.start()

    geo.export_csg(
        title="Converted with GEOUNED",
        geometryName="csg",
        outFormat=(
            "openmc_xml",
            "openmc_py",
            "serpent",
            "phits",
            "mcnp",
        ),
        volSDEF=True,
        volCARD=False,
        UCARD=None,
        dummyMat=True,
        cellCommentFile=False,
        cellSummaryFile=False,
    )

You can also load up JSON configuration files from the API.
The JSON file format is also usable with the Geouned Command Line Interface  and there are more `complete examples of JSON files <python_cli_usage.html>`_ on that section of the documentation.
For this example we assume you have a JSON file called 'config.json' in the same directory as the script.

.. code-block:: python

    import geouned

    geo = geouned.CadToCsg.from_json("config.json")
    geo.start()
    geo.export_csg()
