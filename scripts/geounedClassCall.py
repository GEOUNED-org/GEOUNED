#!/usr/bin/python

# Path to GEOUNED Package


# only if modules are not in the PYTHONPATH or directly installed in the python distribution used
import sys

geo_path = "C:\\Users\\Patrick\\Documents\\GitHub\\GEOUNED\\src"
sys.path.append(geo_path)
sys.path.append("C:\\Program Files\\FreeCAD 0.19\\bin...")

from geouned import CadToCsg

stepFileName = "placa.stp"

GEO = CadToCsg("Conversion Example")

GEO.set("stepFile", stepFileName)
GEO.set("geometryName", "Placa")
GEO.set("outFormat", ("mcnp", "openMC_XML"))
GEO.set("planeDistance", 0.05)
GEO.set("quadricPY", True)
GEO.set("P_abc", "12f")
GEO.set("P_d", "12f")
GEO.Start()
