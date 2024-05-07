Install
=======

There are several options for installing GEOUNED package.
The installation selected has implications for how you run GEOUNED Python scripts.
Currently the Mamba / Conda install is the recommended method.
The main difficulty in installing GEOUNED is integration between the main dependnacy (FreeCAD) and the users system Python.
Mamba / Conda provides the a connection between the FreeCAD Python library and your system Python.
Users have had success installing FreeCAD and making use of the FreeCAD Python as well so these installation methods are listed for completeness.


Linux
-----

Mamba (recommended)
~~~~~~~~~~~~~~~~~~~

First we need to install a Mamba distribution. There are a few options but here we opt for Miniforge3 as it includes Mamba.

You can follow the install instructions for `Miniforge3 here <https://github.com/conda-forge/miniforge>`_ or follows the commands below.
Download 

.. code-block:: sh

    curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh

Install Miniforge3

.. code-block:: sh

    bash Miniforge3-Linux-x86_64.sh


Activate the base environment in your current terminal

.. code-block:: sh

    mamba activate


It is recommended to create a new environment

.. code-block:: sh

    mamba create --name new_new python=3.11


Activate the new environment

.. code-block:: sh

    mamba activate new_new

We have aspirations to create a conda-forge package which will combine these final two steps, but for now FreeCAD and GEOUNED can be installed in two commands.
Install FreeCAD which is the main dependency

.. code-block:: sh

    mamba install -c conda-forge freecad


Install GEOUNED with pip

.. code-block:: sh

    pip install git+https://github.com/GEOUNED-code/GEOUNED.git

Then you will be able to run import GEOUNED from within Python

.. code-block:: python

    import geouned

First we need to install a Conda distribution. There are a few options but we here we opt for `MiniConda3<https://docs.anaconda.com/free/miniconda/>`_ as it downloads quicker than the fuller `AnaConda<https://www.anaconda.com/download>`_.

You can follow the install instructions for `MiniConda3<https://docs.anaconda.com/free/miniconda/>`_ or follow the commands below.
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

    conda create --name new_new python=3.11


Activate the new environment

.. code-block:: sh

    conda activate new_new

We have aspirations to create a conda-forge package which will combine these final two steps, but for now FreeCAD and GEOUNED can be installed in two commands.
Install FreeCAD which is the main dependency

.. code-block:: sh

    conda install -c conda-forge freecad


Install GEOUNED with pip

.. code-block:: sh

    pip install git+https://github.com/GEOUNED-code/GEOUNED.git

Then you will be able to run import GEOUNED from within Python

.. code-block:: python

    import geouned

Apt-get
~~~~~~~

Snap
~~~~

AppImage
~~~~~~~~

Mac
---


Mamba
~~~~~

Conda
~~~~~

Brew
~~~~


Windows
-------

Mamba
~~~~~

Conda
~~~~~

Portable FreeCAD installer
~~~~~~~~~~~~~~~~~~~~~~~~~~

Windows Subsystem for Linux (WSL)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

