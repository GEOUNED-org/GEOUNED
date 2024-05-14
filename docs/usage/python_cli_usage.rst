Command Line Usage
==================

GEOUNED can be used in the command line.

These examples assumes you have a CAD STEP file in the current working directory of the terminal called "cuboid.stp"

The most minimal use case below shows a minimal config.json file being used.

First create a JSON file called "config.json" containing the following.

.. code-block:: json

    {
        "step_file": "cuboid.stp"
    }

Then execute the command line interface tool to convert your sTEP file with the default configuration.

.. code-block:: bash

    geouned_cadtocsg --config config.json

The following example shows a usage with every attributes specified in the config.json file.

#TODO