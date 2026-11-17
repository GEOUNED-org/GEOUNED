Command Line Tool Usage
=======================

GEOUNED CAD to CSG conversion can be performed in the command line.

These examples assumes you have a CAD STEP file in the current working directory of the terminal called "cuboid.stp"

The most minimal use case below shows a minimal config.json file being used.

First create a JSON file called "config.json" containing the following.

.. code-block:: json

    {
        "load_step_file": {
            "filename":"cuboid.stp"
        }
    }

Then execute the command line interface tool to convert your STEP file to CSG files with the default configuration.

.. code-block:: bash

    geouned_cadtocsg -i config.json

The following example shows a usage with every attributes specified in the config.json file.

The contents of the JSON file closely matches the Class arguments and method arguments when using the Python package.

For a full description of each keyword see the `Python API reference section <../python_api.html>`_ of the documentation.

Here is a complete JSON file specification

.. code-block:: json

    {
        "load_step_file": {
            "filename": "cuboid.stp",
            "skip_solids": []
        },
        "Options": {
            "forceCylinder": false,
            "newSplitPlane": true,
            "delLastNumber": false,
            "enlargeBox": 2.0,
            "nPlaneReverse": 0,
            "splitTolerance": 0.0,
            "scaleUp": true,
            "quadricPY": false,
            "Facets": false,
            "prnt3PPlane": false,
            "forceNoOverlap": false
        },
        "Tolerances": {
            "relativeTol": false,
            "relativePrecision": 1e-06,
            "value": 1e-06,
            "distance": 0.0001,
            "angle": 0.0001,
            "pln_distance": 0.0001,
            "pln_angle": 0.0001,
            "cyl_distance": 0.0001,
            "cyl_angle": 0.0001,
            "sph_distance": 0.0001,
            "kne_distance": 0.0001,
            "kne_angle": 0.0001,
            "tor_distance": 0.0001,
            "tor_angle": 0.0001,
            "min_area": 0.01
        },
        "NumericFormat": {
            "P_abc": "14.7e",
            "P_d": "14.7e",
            "P_xyz": "14.7e",
            "S_r": "14.7e",
            "S_xyz": "14.7e",
            "C_r": "12f",
            "C_xyz": "12f",
            "K_xyz": "13.6e",
            "K_tan2": "12f",
            "T_r": "14.7e",
            "T_xyz": "14.7e",
            "GQ_1to6": "18.15f",
            "GQ_7to9": "18.15f",
            "GQ_10": "18.15f"
        },
        "Settings": {
            "matFile": "",
            "voidGen": true,
            "debug": false,
            "compSolids": false,
            "simplify": "no",
            "exportSolids": "",
            "minVoidSize": 200.0,
            "maxSurf": 50,
            "maxBracket": 30,
            "voidMat": [],
            "voidExclude": [],
            "startCell": 1,
            "startSurf": 1,
            "sort_enclosure": false
        },
        "export_csg":{
            "title": "Converted with GEOUNED",
            "geometryName": "csg",
            "outFormat": ["openmc_xml", "openmc_py", "serpent", "phits", "mcnp"],
            "volSDEF": false,
            "volCARD": true,
            "UCARD": null,
            "dummyMat": false,
            "cellCommentFile": false,
            "cellSummaryFile": true
        }
    }

Note that JSON requires ```null``` to be passed in which gets translated to ```None``` in Python.
This is converted in the same way as the minimal JSON config file

.. code-block:: bash

    geouned_cadtocsg -i config.json
