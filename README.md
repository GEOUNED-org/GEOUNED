# GEOUNED
A tool to convert from CAD to CSG & CSG to CAD for Monte Carlo transport codes (MCNP & OpenMC).

## Installation 

Install directly from the repository (for python versions > 3.7 => FreeCAD >= 0.19 ):

```bash
pip install git+https://github.com/GEOUNED-code/GEOUNED.git
```

Otherwise, the source code included in the ~src/ folder can be directly downloaded and properly reached by path variable using as follows:

```python
import sys
GEO_path='the path in your local computer'
sys.path.append('GEO_path')
``` 

the same should be made for FreeCAD libraries. 

If you are using FreeCAD in windows there is included a python distribution within FreeCAD (by default located in C:\Program Files\FreeCAD 0.XX\bin\).
In that case you can install directly de module using:

```bash
C:\Program Files\FreeCAD 0.XX\bin\python.exe -m pip install git+https://github.com/GEOUNED-code/GEOUNED.git
```
with that using this python you have directly access to both FreeCAD and GEOUNED python modules.

## How to use

In the ~scripts/ folder you can find a script (geouned for linux or geouned.py for windows) to call both modules:
* GEOUNED to convert from CAD to CSG for both MCNP and OpenMC codes.
* GEOReverse to convert from CSG to CAD for MCNP.

Detailed descriptions of all options for the user are given in the manual located in ~docs/ folder