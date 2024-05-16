#!/usr/bin/python

# Path to GEOUNED Package


# only if modules are not in the PYTHONPATH or directly installed in the python distribution used
import sys

from geouned import CadToCsg
from geouned.GEOReverse import reverse

# geo_path="C:\\Users\\Juan\\Documents\\work\\GEOUNED\\RepoGit\\GitHub\\GEOUNEDcode\\src"
# sys.path.append(geo_path)
# sys.path.append('C:\\Program Files\\FreeCAD 0.19\\bin...')


runReverse = False
if len(sys.argv) < 2:
    inifile = "config.ini"

elif len(sys.argv) == 2:
    if sys.argv[1] == "-r":
        runReverse = True
        inifile = "configReverse.ini"
    else:
        inifile = sys.argv[1]

elif len(sys.argv) == 3:
    if sys.argv[1] == "-r":
        runReverse = True
        inifile = sys.argv[2]
    elif sys.argv[2] == "-r":
        runReverse = True
        inifile = sys.argv[1]
    else:
        raise ValueError("Bad option")

else:
    raise ValueError("Too many input arguments")

if not runReverse:
    GEO = CadToCsg()
    GEO.set_configuration(inifile)
    GEO.Start()

else:
    print(inifile)
    reverse(inifile)
