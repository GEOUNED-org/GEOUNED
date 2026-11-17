Developers guide
===============

Developer install with Mamba
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First we need to install a Mamba distribution. There are a few options but here we opt for Miniforge3 as it includes Mamba.

Conda could be used instead of Mamba but Mamba is faster and more reliable. If you don't have permission to install Miniforge3 then you could try using Conda instead. For you could install mamba into conda with.

.. code-block:: sh

    conda install -c conda-forge mamba -y

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

As our main dependency FreeCAD is not available on PYPi but we can install it from conda-forge.

.. code-block:: sh

    mamba install -c conda-forge freecad -y

To be able to run test cases and stochastic volume calculation, it is needed to install ``openmc`` and ``pytest`` from conda-forge.

.. code-block:: sh

    mamba install -c conda-forge pytest openmc

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

Keeping your fork tags up to date
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that your fork of the repository should have releases tags which conform to PEP440 i.e. semantic versioning.
When installing, setuptools uses these tags to get the version of GEOUNED.
Release versions in the central repository conform to this standard.
You can keep the tags in your fork up to date with the main repository using:

.. code-block:: sh

   git fetch --tags https://github.com/GEOUNED-org/GEOUNED


Building the docs locally
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: sh
    
        python -m pip install -e .[docs]
        sphinx-build docs _build

Then view the docs by opening the _build/index.html file in a web browser.

When the CI builds docs it puts the latest stable version in the _build directory on the gh-pages branch.

Versions (including dev) are built and put in subdirectories of the _build directory on the gh-pages branch.

Running the tests locally
~~~~~~~~~~~~~~~~~~~~~~~~~

As we installed the tests dependencies using the [tests] option the we can run the tests locally with pytest.


.. code-block:: sh

    python -m pip install -e .[tests,docs]

However we need one more dependency to run the tests.

.. code-block:: sh

    mamba install -c conda-forge openmc -y

Then we can run the tests with the following command from the root of the repository.

.. code-block:: sh

    python -m pytest

We can run individual test files by specifying the file path

.. code-block:: sh

    python -m pytest tests/test_convert.py

We can run individual test functions by specifying the file path and function name

.. code-block:: sh

    python -m pytest tests/test_convert.py -k 'test_conversion'

Additional pytest options that might be useful are including -s for standard output and -vv for very verbose output.

.. code-block:: sh

    python -m pytest -s -vv

Merging a pull requests
~~~~~~~~~~~~~~~~~~~~~~~

Pull requests should be made from feature branches on a fork of the repository to the dev branch.

Tests checking the code will run automatically on the pull request.

If the tests pass and at least one approver approves then the pull request can be merged.

When updating the dev branch from a feature branch then a pull request is should be merged in with the **squashed and merged*** option.

When updating the main branch from the dev branch then the pull request should be merged in with the **create a merge commit** option.

Version numbering
~~~~~~~~~~~~~~~~~

GEOUNED will use Semantic Versioning to number releases of the tool, in the form "Major.Minor.Patch", e.g., “3.15.9”.

Releasing a new version
~~~~~~~~~~~~~~~~~~~~~~~

To release a new version we first need to add and entry to the docs/version_switcher.json file on the dev branch

.. code-block:: python

    [
        {
            "name": "dev",
            "version": "dev",
            "url": "https://geouned-org.github.io/GEOUNED/dev"
        },
        {
            "name": "1.1.0",
            "version": "1.1.0",
            "url": "https://geouned-org.github.io/GEOUNED/1.1.0"
        }
    ]

For example adding version 1.2.3 would look like this

.. code-block:: python

    [
        {
            "name": "dev",
            "version": "dev",
            "url": "https://geouned-org.github.io/GEOUNED/dev"
        },
        {
            "name": "1.1.0",
            "version": "1.1.0",
            "url": "https://geouned-org.github.io/GEOUNED/1.1.0"
        },
        {
            "name": "1.2.3",
            "version": "1.2.3",
            "url": "https://geouned-org.github.io/GEOUNED/1.2.3"
        }
    ]

Then create a `pull request from dev branch to main branch <https://github.com/GEOUNED-org/GEOUNED/compare/main...dev>`_

Once the tests for this pass then merge the pull request in. Use the **create a merge commit** option when merging this pull request from dev to main.

Then `create a new release on the main branch <https://github.com/shimwell/GEOUNED/releases/new>`_ with the version number and a description of the changes.
 
Create a new tag with the version number (e.g. 1.2.3) and the release name (e.g. v1.2.3) and the release description.

Press the Generate release notes button to get the release notes from the pull request descriptions.

Then press the Publish release button to create the release.

This will create the release and trigger github actions for 
- publishing the PyPI package
- building the docs and setting the default docs to the new version

Check the actions both pass by going to the `actions tab https://github.com/shimwell/GEOUNED/actions>`_ on the repository and checking the latest actions.

Conda Forge Releasing
~~~~~~~~~~~~~~~~~~~~~

The conda-forge package release is done after the PyPI release. This is because the conda-forge package is built from the PyPI package.

Conda Forge has a bot that generates a pull request to update the conda-forge recipe. This is done automatically when the PyPI package is released.

The pull request will be generated in the `conda-forge/GEOUNED-feedstock <https://github.com/conda-forge/geouned-feedstock/pulls>`_ repository a short while after the PyPI release.

Check the pull request and if the tests pass then merge the pull request.

A Conda Forge package will be built and released to the conda-forge channel.

Once released the package will be visible on the `conda-forge channel <https://anaconda.org/conda-forge/geouned>`_. 
