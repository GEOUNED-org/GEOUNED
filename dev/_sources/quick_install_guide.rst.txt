Quick Install Guide
===================

Linux
~~~~~

First we need to install a Conda distribution. There are a few options but we here we opt for `MiniConda3 <https://docs.anaconda.com/free/miniconda/>`_ as it downloads quicker than the fuller `Anaconda <https://www.anaconda.com/download>`_.

You can follow the install instructions for `MiniConda3 <https://docs.anaconda.com/free/miniconda/>`_ or follow the commands below.
Download.

.. code-block:: sh

    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh

Install MiniConda3

.. code-block:: sh

    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3


Activate the base environment in your current terminal

.. code-block:: sh

    ~/miniconda3/bin/conda init bash


It is recommended to create a new environment

.. code-block:: sh

    conda create --name new_env python=3.11


Activate the new environment

.. code-block:: sh

    conda activate new_env

Install GEOUNED from conda-forge

.. code-block:: sh

    conda install -c conda-forge geouned -y

Then you will be able to run import GEOUNED from within Python

.. code-block:: python

    import geouned

You will also be able to use the GEOUNED command line tool

.. code-block:: bash

    geouned_cadtocsg --help


Windows
~~~~~~~

First we need to install a Conda distribution. There are a few options but we here we opt for `MiniConda3 <https://docs.anaconda.com/free/miniconda/>`_ as it downloads quicker than the fuller `Anaconda <https://www.anaconda.com/download>`_.

You can follow the install instructions for `MiniConda3 <https://docs.anaconda.com/free/miniconda/>`_ 

Download and execute the Miniforge3 Windows installer

You can get it from the `Miniforge3 GitHub repository here <https://docs.anaconda.com/free/miniconda/>`_ or from the link below

`https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe <https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe>`_

Open "Anaconda Prompt (miniconda3)" which will now be available on your start menu.

It is recommended to create a new environment

.. code-block:: sh

    conda create --name new_env python=3.11

Activate the new environment

.. code-block:: sh

    conda activate new_env

Since one of the GEOUNED dependencies is FreeCAD, classic conda solver of dependencies may significantly slow down the installation. For this reason, it is recommended to install and use conda-libmamba-solver

.. code-block:: sh

    conda install conda-libmamba-solver

Install GEOUNED from conda-forge, forcing the use of libmamba solver

.. code-block:: sh

    conda install -c conda-forge geouned -y --solver=libmamba

Then you will be able to run import GEOUNED from within Python

.. code-block:: python

    import geouned

You will also be able to use the GEOUNED command line tool

.. code-block:: bash

    geouned_cadtocsg --help
