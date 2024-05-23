Linux
=====

User install with Mamba (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    mamba create --name new_env python=3.11


Activate the new environment

.. code-block:: sh

    mamba activate new_env

Install GEOUNED from conda-forge

.. code-block:: sh

    mamba install -c conda-forge geouned

Then you will be able to run import GEOUNED from within Python

.. code-block:: python

    import geouned

You will also be able to use the GEOUNED command line tool

.. code-block:: bash

    geouned_cadtocsg --help

User install with Conda
~~~~~~~~~~~~~~~~~~~~~~~

First we need to install a Conda distribution. There are a few options but we here we opt for `MiniConda3 <https://docs.anaconda.com/free/miniconda/>`_ as it downloads quicker than the fuller `AnaConda <https://www.anaconda.com/download>`_.

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

    conda install -c conda-forge geouned

Then you will be able to run import GEOUNED from within Python

.. code-block:: python

    import geouned

You will also be able to use the GEOUNED command line tool

.. code-block:: bash

    geouned_cadtocsg --help


Developer install with Mamba
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    mamba create --name new_env python=3.11


Activate the new environment

.. code-block:: sh

    mamba activate new_env

We have aspirations to create a conda-forge package which will combine these final two steps, but for now FreeCAD and GEOUNED can be installed in two commands.
Install FreeCAD which is the main dependency

.. code-block:: sh

    mamba install -c conda-forge freecad

Fork the GEOUNED-org/GEOUNED repository by clicking this link, unchecking the Copy the main branch only check box and clicking create fork

`https://github.com/GEOUNED-org/GEOUNED/fork <https://github.com/GEOUNED-org/GEOUNED/fork>`_

Assuming that you have `setup <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent>`_ and `added <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_ SSH keys then we can clone your forked GEOUNED repository.
Replace <USER> with your own github username

.. code-block:: sh

    git clone git@github.com:<USER>/GEOUNED.git

Then change directory into the repository root like this

.. code-block:: sh

    cd GEOUNED

Install GEOUNED with pip, we also prefix this with "python -m" to ensure that pip install uses the correct Python interpreter.
We are also adding the -e to get an editable install so that when you make local changes to the repo these are picked up in your Python scripts straight away (without needing to reinstall).
We also include all the optional dependencies so that we can run tests locally and build the docs locally.

.. code-block:: sh

    python -m pip install -e .[tests,docs]

Then you will be able to run import GEOUNED from within Python

.. code-block:: python

    import geouned

You will also be able to use the GEOUNED command line tool

.. code-block:: bash

    geouned_cadtocsg --help

Checkout feature branches from dev and make local changes on you own branch

.. code-block:: sh

    git checkout dev
    git checkout -b 'my_new_feature'

Pull requests are welcome

.. Apt-get
.. ~~~~~~~

.. Snap
.. ~~~~

.. AppImage
.. ~~~~~~~~

.. Mac
.. ---


.. Mamba
.. ~~~~~

.. Conda
.. ~~~~~

.. Brew
.. ~~~~


.. Windows
.. -------

.. Mamba
.. ~~~~~

.. Conda
.. ~~~~~

.. Portable FreeCAD installer
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~

.. Windows Subsystem for Linux (WSL)
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

