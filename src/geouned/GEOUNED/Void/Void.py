import FreeCAD
import Part

from ..LoadFile import LoadFunctions as LF
from ..Utils.BasicFunctions_part1 import isOposite
from ..Utils.booleanFunction import BoolSequence
from ..Utils.Functions import GeounedSolid, GeounedSurface
from ..Utils.Options.Classes import Options as opt
from ..Void import voidFunctions as VF
from .VoidBoxClass import VoidBox


def voidGeneration(MetaList, EnclosureList, Surfaces, UniverseBox, setting, init):
    voidList = []

    if EnclosureList:
        NestedEnclosure = LF.setEnclosureLevels(EnclosureList)
        VF.assignEnclosure(MetaList, NestedEnclosure)

        # add to Metalist Level 1 enclosures, remove from list material cells totally embedded in Level 1 enclosures
        newMetaList = VF.selectSolids(MetaList, NestedEnclosure[0], UniverseBox)
    else:
        newMetaList = MetaList[:]
        NestedEnclosure = []

    Box = Part.makeBox(
        UniverseBox.XLength,
        UniverseBox.YLength,
        UniverseBox.ZLength,
        FreeCAD.Vector(UniverseBox.XMin, UniverseBox.YMin, UniverseBox.ZMin),
        FreeCAD.Vector(0, 0, 1),
    )

    EnclosureBox = GeounedSolid(None, Box)
    if setting["voidMat"]:
        voidMat = setting["voidMat"]
        EnclosureBox.setMaterial(voidMat[0], voidMat[1], voidMat[2])

    # get voids in 0 Level Enclosure (original Universe)
    # if exist Level 1 enclosures are considered as material cells
    print("Build Void highest enclosure")

    voids = GetVoidDef(newMetaList, Surfaces, EnclosureBox, setting, Lev0=True)
    voidList.append(voids)

    # Perform enclosure void
    # Loop until the lowest enclosure level

    for i, Level in enumerate(NestedEnclosure):

        print("Build Void highest enclosure")
        for j, encl in enumerate(Level):
            if encl.CellType == "envelope":
                continue
            newMetaList = VF.selectSolids(MetaList, encl.SonEnclosures, encl)
            print("Build Void enclosure {} in enclosure level {}".format(j, i + 1))
            # select solids overlapping current enclosure "encl", and lower level enclosures
            voids = GetVoidDef(newMetaList, Surfaces, encl, setting)
            voidList.append(voids)

    voidList.append(setGraveyardCell(Surfaces, UniverseBox))

    return VF.updateVoidList(init, voidList, NestedEnclosure, setting["sortEnclosure"])


def GetVoidDef(MetaList, Surfaces, Enclosure, setting, Lev0=False):

    maxsurf = setting["maxSurf"]
    maxbracket = setting["maxBracket"]
    minSize = setting["minVoidSize"]

    if "full" in setting["simplify"].lower():
        simplifyVoid = "full"
    elif "void" in setting["simplify"].lower():
        simplifyVoid = "diag"
    else:
        simplifyVoid = "no"

    if Lev0:
        Universe = VoidBox(MetaList, Enclosure.BoundBox)
    else:
        Universe = VoidBox(MetaList, Enclosure.CADSolid)

    Initial = [Universe]
    VoidDef = []
    iloop = 0
    while iloop < 50:
        Temp = []
        iloop += 1
        nvoid = len(Initial)
        print("Loop, Box to Split :", iloop, nvoid)

        for iz, z in enumerate(Initial):
            nsurfaces, nbrackets = z.getNumbers()
            if opt.verbose:
                print(
                    "{} {}/{} {} {}".format(iloop, iz + 1, nvoid, nsurfaces, nbrackets)
                )

            if nsurfaces > maxsurf and nbrackets > maxbracket:
                newspace = z.Split(minSize)
            else:
                newspace = None

            if type(newspace) is tuple:
                Temp.extend(newspace)
            else:
                #           if len(z.Objects) >= 50 : z.refine()
                boxDim = (
                    z.BoundBox.XMin * 0.1,
                    z.BoundBox.XMax * 0.1,
                    z.BoundBox.YMin * 0.1,
                    z.BoundBox.YMax * 0.1,
                    z.BoundBox.ZMin * 0.1,
                    z.BoundBox.ZMax * 0.1,
                )

                print("build complementary {} {}".format(iloop, iz))

                cell, CellIn = z.getVoidComplementary(Surfaces, simplify=simplifyVoid)
                if cell is not None:
                    VoidCell = (cell, (boxDim, CellIn))
                    VoidDef.append(VoidCell)

        Initial = Temp
        if len(Temp) == 0:
            break

    voidList = []
    for i, vcell in enumerate(VoidDef):
        mVoid = GeounedSolid(i)
        mVoid.Void = True
        mVoid.CellType = "void"
        mVoid.setDefinition(vcell[0], simplify=True)
        mVoid.setMaterial(Enclosure.Material, Enclosure.Rho, Enclosure.MatInfo)
        mVoid.setDilution(Enclosure.Dilution)

        mVoid.__commentInfo__ = vcell[1]

        voidList.append(mVoid)

    return voidList


def setGraveyardCell(Surfaces, UniverseBox):
    Universe = VoidBox([], UniverseBox)

    externalBox = getUniverseComplementary(Universe, Surfaces)
    center = UniverseBox.Center
    radius = 0.51 * UniverseBox.DiagonalLength
    sphere = GeounedSurface(("Sphere", (center, radius)), UniverseBox)
    id, exist = Surfaces.addSphere(sphere)

    sphdef = BoolSequence(str(-id))
    sphdef.operator = "AND"
    sphdef.append(externalBox)

    notsph = BoolSequence(str(id))

    mVoidSphIn = GeounedSolid(0)
    mVoidSphIn.Void = True
    mVoidSphIn.CellType = "void"
    mVoidSphIn.setDefinition(sphdef)
    mVoidSphIn.setMaterial(0, 0, "Graveyard_in")
    mVoidSphIn.__commentInfo__ = None

    mVoidSphOut = GeounedSolid(1)
    mVoidSphOut.Void = True
    mVoidSphOut.CellType = "void"
    mVoidSphOut.setDefinition(notsph)
    mVoidSphOut.setMaterial(0, 0, "Graveyard")
    mVoidSphOut.__commentInfo__ = None

    return (mVoidSphIn, mVoidSphOut)


def getUniverseComplementary(Universe, Surfaces):
    Def = BoolSequence(operator="OR")
    for p in Universe.getBoundPlanes():
        id, exist = Surfaces.addPlane(p)
        if not exist:
            Def.elements.append(-id)
        else:
            s = Surfaces.getSurface(id)
            if isOposite(p.Surf.Axis, s.Surf.Axis):
                Def.elements.append(id)
            else:
                Def.elements.append(-id)
    return Def


def voidCommentLine(CellInfo):
    boxDef, cellIn = CellInfo
    cells = ", ".join(map(str, cellIn))
    box = ", ".join(f"{num:.3f}" for num in boxDef)
    line = f"Automatic Generated Void Cell. Enclosure({box})\n"
    line += f"Enclosed cells : ({cells})\n"
    return line
