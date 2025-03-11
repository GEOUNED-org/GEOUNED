Example of use
==============

The following is an example of how to use GEOUNED to convert a CAD file to a CSG file 
for use in MCNP. It includes the methodology and rationale behind the process.

.. warning:: It is assumed that GEOUNED has been properly installed as described in the 
    `Quick install guide <./quick_install_guide.html>`_.

1. `CAD simplification`_
2. `Export as STEP`_
3. `Conversion to CSG`_
4. `Suspicious solids and debug`_
5. `Void check verification`_

.. _CAD simplification:

CAD simplification
~~~~~~~~~~~~~~~~~~

The CAD model of the system should undergo a process of simplification. In this process,
features that do not affect the radiation transport should be removed. The level of 
detail is usually decreased. Surfaces that are not allowed like splines are substituted
by simpler surfaces like planes or cylinders.

.. _Export as STEP:

Export as STEP
~~~~~~~~~~~~~~

Once the CAD is simplified, it is exported as a STEP file to be used as input by 
GEOUNED.

.. _Conversion to CSG:

Conversion to CSG
~~~~~~~~~~~~~~~~~

The below code is run in a Python environment with GEOUNED installed. In the comments
of the code there is an explanation of the parameters used.

.. code-block:: python

    import geouned

    settings = geouned.Settings(
        debug=True,  # This will create a folder with the CAD bodies and their decomposition
        startCell=12001,  # The first cell id
        startSurf=12001,  # The first surface id
        compSolids=False,  # False so the different solids of a component 
                           # are not combined into a single MCNP cell
        voidGen=True,  # Generate the void cells 
    )

    geo = geouned.CadToCsg(settings=settings)
    geo.load_step_file(filename=f"Geometry/model.stp")  # Path to the STEP file
    geo.start()
    print("Start exporting")
    geo.export_csg(
        geometryName="My model",  # Name of the model to appear in the MCNP file
        outFormat=("mcnp",),  # The output format
        UCARD=50,  # Makes the model a filler universe with filler id 50
        volCARD=True,  # The volume card is included in the cell definition
        volSDEF=True,  # Generate a SDEF card
    )
    print("Exporting finished")

.. _Suspicious solids and debug:

Suspicious solids and debug
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running the above code will generate a debug folder that contains each of the converted 
solids, both the original and the decomposed ones. Those bodies can be opened and 
examined in a CAD software to compare them to the original model. If the decomposed 
model doesn't look right, it is possible to manually decompose a specific solid and 
redo the conversion.

Also, it is possible that a folder called *Suspicious solids* is created. This folder 
contains the same solids as the debug folder but only the ones that have been flagged as
potentially problematic. This can be useful to quickly identify the problematic solids.

After applying any fix necessary and exporting to STEP, the conversion can be run again.

.. _Void check verification:

Void check verification
~~~~~~~~~~~~~~~~~~~~~~~

Once the conversion is finished, it is recommended to stocastically check the volume
of the solids and to make some plots.

In this example, the *UCARD* was used to make the model a filler universe. For these 
initial checks the *U=50* cards should be commented out as we want to test the 
standalone model. Plots with the MCNP plotter should be done to check that the model 
looks as expected. 

Also, a void check simulation should be run. The *volSDEF=True* card automatically 
generated a spherical inward-oriented source with the correct weight to calculate 
volumes stochastically. The tallies for the solid MCNP cells were also automatically 
generated. The user can run the simulation directly and the results of the tallies 
should all be close to 1.0. If a MCNP cell tally shows a result far from 1.0 it would 
mean that the volume of the solid is not correct.

.. warning:: The void check simulation should be run with a high number of particles to 
    achieve an acceptable statistical error. An usually number is 1e9 particles. Always
    check that the statistical error of the tallies is below 0.1.
