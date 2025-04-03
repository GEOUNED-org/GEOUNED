from importlib.metadata import version

# this try except attempts to import freecad (lowercase) which is the conda
# package name for FreeCAD (mixed case) upon import the conda package appends
# the sys path for Conda installed FreeCAD, consequently FreeCAD can then be
# found by subsequent import statements through out the code base
try:
    import freecad
except ImportError:
    pass

from .GEOReverse import *
from .GEOUNED import *

__version__ = version("geouned")

__all__ = ["__version__"]
