[![CI testing](https://github.com/GEOUNED-org/GEOUNED/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/GEOUNED-org/GEOUNED/actions/workflows/ci.yml)
[![Upload Python Package](https://github.com/GEOUNED-org/GEOUNED/actions/workflows/python-publish.yml/badge.svg)](https://github.com/GEOUNED-org/GEOUNED/actions/workflows/python-publish.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)
[![documentation](https://github.com/GEOUNED-org/GEOUNED/actions/workflows/documentation.yml/badge.svg)](https://github.com/GEOUNED-org/GEOUNED/actions/workflows/documentation.yml)

:book: [Documentation](https://GEOUNED-org.github.io/GEOUNED/)


# GEOUNED
A tool to convert from CAD to CSG & CSG to CAD for Monte Carlo transport codes (MCNP & OpenMC).
This repository contains the implementation of the algorithm presented in the paper [GEOUNED: A new conversion tool from CAD to Monte Carlo geometry](https://doi.org/10.1016/j.net.2024.01.052).

## Installation 

Install directly from the repository (for python versions > 3.7 => FreeCAD >= 0.19 ):

```bash
pip install git+https://github.com/GEOUNED-code/GEOUNED.git
```

Otherwise, the source code included in the ~src/ folder can be directly downloaded and properly reached by path variable as follows:

```python
import sys
GEO_path='the path in your local computer'
sys.path.append('GEO_path')
```

the same should be made for FreeCAD libraries. You can also define appropriately the PYTHONPATH variable for both modules. 

If you are using FreeCAD in windows there is included a python distribution within FreeCAD distribution (by default located in C:\Program Files\FreeCAD 0.XX\bin\).
In that case you can install directly de module using:

```bash
C:\Program Files\FreeCAD 0.XX\bin\python.exe -m pip install git+https://github.com/GEOUNED-code/GEOUNED.git
```
using this option you have directly access to both FreeCAD and GEOUNED python modules. 
Furthermore, using this python compatibilities problems between different versions of python are avoided (some dynamic libraries of FreeCAD depends on the version of python used during the built process). 

## How to use

The code is used via python scripting.  

The first step is to call the python modules of both FreeCAD and GEOUNED (this can be avoided if you have installed GEOUNED as commented in the previous section and using the python distributed with FreeCAD).

```python
import sys
GEO_path='the path in your local computer, the complete or relative path to the ~src/ folder'
FreeCAD_path='path of FreeCAD python modules'
sys.path.append(GEO_path)
sys.path.append(FreeCAD_path)
``` 
The second step is to call or GEOUNED for CAD to CSG conversion or GEOReverse for CSG to CAD.

From CAD to CSG:

```python
from geouned import CadToCsg
inifile='Name of config file for forward conversion'
GEO = CadToCsg(inifile)
GEO.set_options()
GEO.start()
``` 
From CSG to CAD (so far only for MCNP): 

```python
from GEOReverse.reverse import reverse
inifile='Name of config file for reverse conversion'
reverse(inifile)
``` 
In the ~scripts/ folder you can find a script (geouned for linux or geouned.py for windows) to call both modules directly by command line as follows:

windows (forward and reverse):
```bash
/~>python geouned.py file_config_name
/~>python geouned.py -r file_config_name
```

linux (forward and reverse):
```bash
/~>geouned file_config_name
/~>geouned -r file_config_name
```

Detailed descriptions of all options of the config files (forward and reverse) are given in the manual located in ~docs/ folder.

## Citation

```
@article{CATALAN2024,
  author = {J.P. Catalán and P. Sauvan and J. García and J. Alguacil and F. Ogando and J. Sanz},
  title = {GEOUNED: A new conversion tool from CAD to Monte Carlo geometry},
  journal = {Nuclear Engineering and Technology},
  year = {2024},
  issn = {1738-5733},
  doi = {https://doi.org/10.1016/j.net.2024.01.052}
}
```
