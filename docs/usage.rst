Usage
=====

GEOUNED can be used as a Python package or via the command line.

Python package usage
--------------------

The Python API has two main classes.
The first main class is ``CadToCsg()`` which converts CAD geometry to Constructive Solid Geometry (CSG).
There are many arguments that can be passed into the ``CadToCsg()`` class which are documented in the Python API section.


The most minimal use case below shows GEOUNED being imported and the CadToCsg being used to convert a STEP CAD file called 'cuboid.stp' into a varity of CSG format. 
If you have install GEOUNED and FreeCAD into your system PYthon then you can simply run a .py script with python.


.. code-block:: python

    import geouned
    geo = geouned.CadToCsg(step_file = 'cuboid.stp')
    geo.start()
    geo.export_csg()

The above examples makes use of default :meth:`geouned.Options`, :meth:`geouned.Tolerances`, :meth:`geouned.Settings`, :meth:`geouned.NumericFormat` and :meth:`geouned.Settings`.
Users can change any of these to suit the conversion desired.
The following example changes several default values of the conversion.
The example also changes some of the defaults in the :meth:`geouned.CadToCsg.export_csg` to control how the CSG file is written.

.. code-block:: python

    import geouned

    my_options = geouned.Options(
        forceCylinder=True,
        splitTolerance=0,
        newSplitPlane=True,
        nPlaneReverse=0,
    )

    my_settings = geouned.Settings(
        compSolids=False,
        voidGen=True,
        minVoidSize=100,
        debug=False,
        simplify="no",
    )

    geo = geouned.CadToCsg(
        stepFile=f"{input_step_file.resolve()}", settings=my_settings, options=my_options
    )

    geo.start()

    geo.export_csg(
        title="Input Test",
        geometry_name="test",
        out_formats=("mcnp", "openmc_xml", "openmc_py", "serpent", "phits"),
        volCARD=False,
        volSDEF=True,
        cellSummaryFile=False,
        cellCommentFile=False,
        dummyMat=True,
    )

GEOUNED can also accept JSON or `JSON5 <https://json5.org/>`_ format input files that can contain all the configuration options.
A file in JSON or `JSON5 <https://json5.org/>`_ format can be read in and used to set the CadToCsg instance attributes and the CadToCsg.export_csg arguments.
Here is an example `JSON5 <https://json5.org/>`_ format file.

.. code-block:: python

    {
        step_file: "stepfilename.stp",
        Settings:{
            matFile: "materials.txt",
            compSolids: true, 
            startCell:1,
        },
        export_csg:{
            title : "title of the model in MCNP input",
            UCARD : 101,
            outFormat: ["mcnp", "openmc_py", "openmc_xml"],
            geometryName: "pieza",
            volSDEF  :false,
            volCARD  :false,
            dummyMat: false,
        }
    }

Here is the Python code to convert the cad using this config file 

.. code-block:: python

    import geouned

    # this loads the config which contains all the class args and the args needed
    # for an export call and does the whole conversion
    geouned.CadToCsg.from_config("config.json5")



Command line usage
------------------

GEOUNED can also be used from the command line with the command

This example loads up a config file with the default name (config.json)

.. code-block:: bash

    geouned_cadtocsg

This example specifies a JSON config file name called "example.json"

.. code-block:: bash

    geouned_cadtocsg -i configuration.json

Both JSON and `JSON5 <https://json5.org/>`_ file formats are supported for config files.
This example specifies a JSON5 config file name called "example.json5"

.. code-block:: bash

    geouned_cadtocsg -i configuration.json5

This example prints the options to the command line.

.. code-block:: bash

    geouned_cadtocsg --help
