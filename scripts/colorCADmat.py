import re
import sys

sys.path.append("/usr/lib64/freecad/lib64/")
import ImportGui

bd = 255
gd = 255 * 255
rd = 255 * 255 * 255
matNumber = re.compile(r"Material_(?P<matnum>\d+)")


def getColor(num, vmin=0, vmax=rd):
    colRange = rd
    dnum = vmax - vmin
    dx = num - vmin
    scale = colRange / dnum
    snum = dx * scale + 1

    red = snum // gd
    r = snum - red * gd
    green = r // bd
    blue = r - green * bd
    return (red / 255, green / 255, blue / 255)


def setColorMaterial(documents):
    featureObj = []
    matMin = 999999999
    matMax = 0
    for obj in documents.Objects:
        if obj.TypeId == "Part::Feature":
            for o in obj.InListRecursive:
                m = matNumber.search(o.Label)
                if m:
                    mat = int(m.group("matnum"))
                    matMin = min(matMin, mat)
                    matMax = max(matMax, mat)
                    featureObj.append((obj, mat))

    for obj, mat in featureObj:
        obj.ViewObject.ShapeColor = getColor(mat, matMin, matMax)


# color solids in a Freecad document with respect to the material number
doc = App.ActiveDocument
setColorMaterial(doc)
name = doc.Label
outname = name + "_color"
ImportGui.export(doc.RootObjects, outname + ".stp")
doc.saveAs(outname + ".FCStd")
