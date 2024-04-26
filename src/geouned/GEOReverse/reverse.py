import FreeCAD
import Import

from .CodeVersion import *
from .Modules.buildCAD import buildCAD, makeTree
from .Modules.McnpInput import McnpInput
from .Modules.Objects import CADCell
from .Modules.processInp import setOptions
from .Modules.XmlInput import XmlInput


def reverse(optFile="configRevese.ini"):

    printCodeVersion()

    setting = setOptions(optFile)

    geomfile = setting["fileIn"]
    outname = setting["fileOut"]
    outBox = setting["outBox"]
    inFormat = setting["inFormat"]

    CADselection = {
        "Ustart": setting["UStart"],
        "levelMax": setting["levelMax"],
        "cell": setting["cell"],
        "mat": setting["mat"],
        "format": setting["inFormat"],
    }

    UnivCell = CadCell()
    UnivCell.shape = UnivCell.makeBox(FreeCAD.BoundBox(*outBox))

    # get geometry definition from MCNP input
    if inFormat == "mcnp":
        geom = McnpInput(geomfile)
    elif inFormat == "openMC_XML":
        geom = XmlInput(geomfile)
    else:
        msg = (
            f"input format type {inFormat} is not supported."
            'Supported options are "mcnp" or "openMC_XML"'
        )
        raise ValueError(msg)

    CADCells, fails = buildCAD(UnivCell, geom, CADselection)

    if fails:
        print("failed in conversion", fails)

    CADdoc = FreeCAD.newDocument("WorkingDoc")

    makeTree(CADdoc, CADCells)
    Import.export(CADdoc.Objects[0:1], outname + ".stp")
    CADdoc.saveAs(outname + ".FCStd")


def printCodeVersion():

    FreeCAD_Version = "{V[0]:}.{V[1]:}.{V[2]:}".format(V=FreeCAD.Version())
    title = """\
#########################################################################
#                                                                       # 
#      GEOReverse version {:<11}{}{:>26}
#      FreeCAD    version {:<11}{:>36}  
#                                                                       # 
#########################################################################""".format(
        GEOReverse_Version, GEOReverse_ReleaseDate, "#", FreeCAD_Version, "#"
    )
    print(title)


if __name__ == "__main__":
    reverse()
