Command Line Tool Usage
=======================

GEOUNED CSG to CAD conversion can be performed in the command line.

Both OpenMC XML CSG and MCNP CSG formats are supported.

The first example assumes you have an OpenMC XML CSG file in the current working directory of the terminal called "cylinder_box.xml".

The most minimal use case below shows a minimal config_openmc.json file being used.

First create a JSON file called "config_openmc.json" containing the following.

.. code-block:: json

    {
        "export_cad":{
            "input_filename": "cylinder_box.xml",
            "csg_format": "openmc_xml"
        }
    }


Then execute the command line interface tool to convert your OpenMC XML CSG file to a STEP CAD file with the default configuration.

.. code-block:: bash

    geouned_csgtocad -i config_openmc.json

MCNP CSG files can also be converted, this example assumes you have a MCNP CSG file in te current working directory called "cylinder_box.mcnp".

For MCNP CSG files, the JSON config file should be as follows and is assumed to be named "config_mcnp.json".


.. code-block:: json

    {
        "export_cad":{
            "input_filename": "cylinder_box.mcnp",
            "csg_format": "mcnp"
        }
    }

Then execute the command line interface tool to convert your MCNP CSG file to a STEP CAD file with the default configuration.

.. code-block:: bash

    geouned_csgtocad -i config_mcnp.json


The following example shows a usage with every attributes specified in the config.json file.

The contents of the JSON file closely matches the Class arguments and method arguments when using the Python package.

For a full description of each keyword see the `Python API reference section <../python_api.html>`_ of the documentation.

Here is a complete JSON file specification

.. code-block:: json

    {
    "export_cad":{
        "input_filename": "tests/csg_files/cylinder_box.xml",
        "csg_format": "openmc_xml",
        "output_filename" : "cad_from_csg",
        "bounding_box": [-1000, -1000, -1000, 1000, 1000, 1000],
        "universe_start": 0,
        "level_max": "all",
        "cell_range_type": "all",
        "cell_range": [],
        "mat_range_type": "all",
        "mat_range": []
        }
    }

.. code-block:: bash

    geouned_csgtocad -i config.json
