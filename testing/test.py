import os
import subprocess
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parents[1] / "src"))
sys.path.append("/usr/lib64/freecad/lib64/")

from geouned import CadToCsg


def setInput(inName, inpDir, outDir):

    if inName.endswith(".step"):
        filename = inName[0:-5]
    elif inName.endswith(".stp"):
        filename = inName[0:-4]
    else:
        filename = inName

    if inpDir == "":
        inpDir = "."
    if outDir == "":
        outDir = "."

    inName = f"{inpDir}/{inName}"
    outName = f"{outDir}/{filename}"

    template = """[Files]
title    = Input Test
stepFile = {}
geometryName = {}

[Parameters]
compSolids = False
volCARD    = False
volSDEF    = True
voidGen    = True
dummyMat    = True
minVoidSize =  100
cellSummaryFile = False
cellCommentFile = False
debug       = False
simplify   = full

[Options]
forceCylinder = False
splitTolerance = 0
newSplitPlane = True
nPlaneReverse = 0
""".format(
        inName, outName
    )

    with open(file="config.ini", mode="w", encoding="utf-8") as outfile:
        outfile.write(template)


def getInputList(folder, ext=None):
    filenameList = []
    if ext is None:
        return os.listdir(folder)

    elif type(ext) is str:
        if ext[0] != ".":
            ext = "." + ext
        for f in os.listdir(folder):
            if f.endswith(ext):
                filenameList.append(f)

    elif type(ext) is list or type(ext) is tuple:
        filenameList = []
        extList = []
        for e in ext:
            if e[0] != ".":
                e = "." + e
            extList.append(e)

        for f in os.listdir(folder):
            for e in extList:
                if f.endswith(e):
                    filenameList.append(f)

    return filenameList


def runMCNP(path, inpFile):
    code = "/opt/mcnp5/last/bin/mcnp5"
    xsdir = "/opt/mcnp.data/ascii/xsdir"
    inp = inpFile
    out = inpFile[0:-1] + "o"
    mctal = inpFile[0:-1] + "m"
    cmd = "cd {} && {} i={} o={} mctal={} xsdir={}".format(path, code, inp, out, mctal, xsdir)
    os.system(cmd)


def clean(folder):
    filename = "{}/runtpe".format(folder)
    if os.path.isfile(filename):
        os.remove(filename)


def cleanDir(folder):
    for file in getInputList(folder, (".o", ".m")):
        filename = "{}/{}".format(folder, file)
        os.remove(filename)


def checkLost(outp):
    cmd = "grep 'particles got lost.' {}".format(outp)
    stdout = subprocess.getoutput(cmd)
    lineout = stdout.split("\n")
    lost = 0
    for l in lineout:
        if l == "" or "terminated" in l:
            continue
        lost = int(l[0 : l.index("particles")])
    return lost


def getMctalValues(mctal):
    file = open(mctal, "r")
    while True:
        line = file.readline()
        if "vals" in line:
            break

    rawData = []
    line = file.readline()
    while "tfc" not in line:
        rawData.extend(line.split())
        line = file.readline()
    file.close()

    values = []
    for i in range(int(len(rawData) / 2)):
        val = float(rawData[2 * i])
        err = float(rawData[2 * i + 1])
        values.append((val, err))
    return values


def printResults(f, res, lost):

    line = "{} :\n".format(f)
    if lost != 0:
        line += "{} particles got lost.\n".format(lost)

    for i, vals in enumerate(res):
        line += "  {:<3d} : {:8.5f}  {:7.5f}\n".format(i + 1, vals[0], vals[1])
    print(line)


def mkGEOInp(inpDir, outDir):
    for f in getInputList(inpDir, ("stp", "step")):
        setInput(f, inpDir, outDir)
        GEO = CadToCsg()
        GEO.set_configuration(inifile)
        GEO.Start()
        del GEO


def processMCNPfolder(outDir):
    cleanDir(outDir)
    for f in getInputList(outDir, ".mcnp"):
        os.rename(f"{outDir}/{f}", f"{outDir}/{f[:-4]}i")

    for f in getInputList(outDir, ".i"):
        runMCNP(outDir, f)
        clean(outDir)


def postProcess(folder):
    for fmct in getInputList(folder, ".m"):
        fout = fmct[:-1] + "o"
        mctal = "{}/{}".format(folder, fmct)
        outp = "{}/{}".format(folder, fout)
        res = getMctalValues(mctal)
        lost = checkLost(outp)
        printResults(fmct[:-1] + "i", res, lost)


# ****************************************************

misInp = "inputSTEP/Misc"
misOut = "outMCNP/Misc"
cylInp = "inputSTEP/DoubleCylinder"
cylOut = "outMCNP/DoubleCylinder"
torInp = "inputSTEP/Torus"
torOut = "outMCNP/Torus"
lrgInp = "inputSTEP/large"
lrgOut = "outMCNP/large"
ms0Inp = "inputSTEP/"
ms0Out = "outMCNP/"

inifile = "config.ini"
folder = {
    "misc": (misInp, misOut),
    "cyl": (cylInp, cylOut),
    "torus": (torInp, torOut),
    "large": (lrgInp, lrgOut),
    "misc0": (ms0Inp, ms0Out),
}

test = "all"

if test == "all":
    for inpDir, outDir in folder.values():
        mkGEOInp(inpDir, outDir)

    for inpDir, outDir in folder.values():
        processMCNPfolder(outDir)

    for inpDir, outDir in folder.values():
        postProcess(outDir)

else:
    inpDir, outDir = folder[test]
    mkGEOInp(inpDir, outDir)
    processMCNPfolder(outDir)
    postProcess(outDir)
