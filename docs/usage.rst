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
    geo = geouned.CadToCsg(
        step_file = 'cuboid.stp',
        out_formats = ['openmc_xml', 'openmc_py', 'PHITS', 'Serpent', 'MCNP']
    )
    geo.start()


Command line usage
------------------

TODO