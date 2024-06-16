.. usage

Install
=======

There are several options for installing GEOUNED package.
The installation selected has implications for how you run GEOUNED Python scripts.
Currently the Mamba / Conda install is the recommended method.
The main complication when installing GEOUNED is integration between the main dependency (FreeCAD) and the users system Python.
Mamba / Conda provides a connection between the FreeCAD Python library and your system Python.
Users have also had success installing FreeCAD and making use of the Python version inbuilt into FreeCAD and the freecad.cmd however the integration of FrreCAD with the system Python is more challenging when installing in this manner.
.. TODO as well so these installation methods are listed for completeness.

Note that your fork of the repostory must have releases which conform to PEP440 i.e. semantic versioning. When installing, setuptools
uses these tags to get the release number of GEOUNED. Release versions in the central repository now conform to this standard and can be reflected in your
fork using:

.. code-block:: sh

   git fetch --tags upstream

where upstream is the remote url of the original repository.

.. toctree::
   :numbered:
   :maxdepth: 1

   install_linux
   install_windows