Python Package Usage, CSG to CAD conversion
===========================================

The main class is ``CsgToCad()`` which converts Constructive Solid Geometry (CSG) to CAD.
There are a few arguments that can be passed into the ``CsgToCad()`` class which are documented in the `Python API reference section <../python_api.html>`_ of the documentation.


If you have install GEOUNED and FreeCAD into your system Python then you can simply run a .py script with Python.
The most minimal use case below shows GEOUNED being imported and the CsgToCad being used to convert a CSG geometry into a STEP CAD file. 
The example makes use of default  attributes.

.. code-block:: python

    import geouned


Users can change :meth:`geouned.Options`, :meth:`geouned.Settings`, :meth:`geouned.Tolerances` and :meth:`geouned.NumericFormat` to suit the conversion desired.
The following example shows a usage with every attributes specified.

.. code-block:: python

    import geouned
