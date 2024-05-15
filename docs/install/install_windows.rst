Windows
=======


User install with Mamba (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First we need to install a Mamba distribution. There are a few options but here we opt for Miniforge3 as it includes Mamba.

You can follow the install instructions for `Miniforge3 here <https://github.com/conda-forge/miniforge>`_ or follows the commands below.

Download and execute the Miniforge3 Windows installer

You can get it from the Miniforge3 GitHub repository `here <https://github.com/conda-forge/miniforge?tab=readme-ov-file#miniforge-pypy3>`_ or from the link below

`https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>`_

Follow the prompts and complete the installation process

Open "Miniforge Prompt" which will now be available on your start menu.


It is recommended to create a new environment

.. code-block:: sh

    mamba create --name new_env python=3.11

Activate the new environment

.. code-block:: sh

    mamba activate new_env

We have aspirations to create a conda-forge package which will combine these final two steps, but for now FreeCAD and GEOUNED can be installed in two commands.
Install FreeCAD which is the main dependency

.. code-block:: sh

    mamba install -c conda-forge freecad


Install GEOUNED with pip, we also prefix this with "python -m" to ensure that pip install uses the correct Python interpreter.

.. code-block:: sh

    python -m pip install geouned

Then you will be able to run import GEOUNED from within Python

.. code-block:: python

    import geouned

You will also be able to use the GEOUNED command line tool

.. code-block:: bash

    geouned_cadtocsg --help
