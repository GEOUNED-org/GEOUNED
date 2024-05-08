import logging

# this try except attempts to import freecad (lowercase) which is the conda
# package name for FreeCAD (mixed case) upon import the conda package appends
# the sys path for Conda installed FreeCAD, consequently FreeCAD can then be
# found by subsequent import statements through out the code base
try:
    import freecad
except ImportError:
    pass

from .GEOUNED import *
from .GEOReverse import *

logging.basicConfig(
    filename="geouned.log",
    filemode="w",
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s :: %(levelname)s :: %(funcName)s :: %(lineno)d :: %(message)s",
)
