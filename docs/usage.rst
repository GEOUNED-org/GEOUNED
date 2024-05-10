Usage
=====

GEOUNED can be used as a Python package or via the command line.

Python package usage
--------------------

The Python API has two main classes.
The first main class is ``CadToCsg()`` which converts CAD geometry to Constructive Solid Geometry (CSG).
There are many arguments that can be passed into the ``CadToCsg()`` class which are documented in the `Python API section <python_api.html>`_ of the documentation.


If you have install GEOUNED and FreeCAD into your system Python then you can simply run a .py script with Python.
The most minimal use case below shows GEOUNED being imported and the CadToCsg being used to convert a STEP CAD file called 'cuboid.stp' into a vanity of CSG format. 
The example makes use of default  attributes.

.. code-block:: python

    import geouned
    geo = geouned.CadToCsg(stepFile = 'cuboid.stp')
    geo.start()
    geo.export_csg()

Users can change :meth:`geouned.Options`, :meth:`geouned.Tolerances` and :meth:`geouned.NumericFormat` to suit the conversion desired.
The following example shows a more complete usage with several default attributes changed.

.. code-block:: python

    import geouned

    my_options = geouned.Options(
        forceCylinder=True,
        splitTolerance=0,
        newSplitPlane=True,
        nPlaneReverse=0,
    )

    my_tolerances = geouned.Tolerances(
        min_area=0.011
    )
    my_numeric_format = geouned.NumericFormat(
        C_r="13f"
    )

    geo = geouned.CadToCsg(
        stepFile='cuboid.stp',
        options=my_options,
        tolerances=my_tolerances
        numeric_format=my_numeric_format
    )

    geo.start()
