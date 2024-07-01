import logging

import FreeCAD
import Part
from tqdm import tqdm

from ..loadfile import load_functions as LF
from ..utils.basic_functions_part1 import is_opposite
from ..utils.boolean_function import BoolSequence
from ..utils.functions import GeounedSolid, GeounedSurface
from . import void_functions as VF
from .void_box_class import VoidBox

logger = logging.getLogger("general_logger")


def void_generation(
    MetaList,
    EnclosureList,
    Surfaces,
    UniverseBox,
    setting,
    init,
    options,
    tolerances,
    numeric_format,
):
    voidList = []

    if EnclosureList:
        NestedEnclosure = LF.set_enclosure_levels(EnclosureList)
        VF.assignEnclosure(MetaList, NestedEnclosure)

        # add to Metalist Level 1 enclosures, remove from list material cells totally embedded in Level 1 enclosures
        newMetaList = VF.select_solids(MetaList, NestedEnclosure[0], UniverseBox)
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
    if setting.voidMat:
        voidMat = setting.voidMat
        EnclosureBox.set_material(voidMat[0], voidMat[1], voidMat[2])

    # get voids in 0 Level Enclosure (original Universe)
    # if exist Level 1 enclosures are considered as material cells
    logger.info("Build Void highest enclosure")

    voids = get_void_def(
        newMetaList,
        Surfaces,
        EnclosureBox,
        setting,
        options,
        tolerances,
        numeric_format,
        Lev0=True,
    )
    voidList.append(voids)

    # Perform enclosure void
    # Loop until the lowest enclosure level

    for i, Level in enumerate(NestedEnclosure):

        logger.info("Build Void highest enclosure")
        for j, encl in enumerate(Level):
            if encl.CellType == "envelope":
                continue
            newMetaList = VF.select_solids(MetaList, encl.SonEnclosures, encl)
            logger.info(f"Build Void enclosure {j} in enclosure level {i + 1}")
            # select solids overlapping current enclosure "encl", and lower level enclosures
            voids = get_void_def(
                newMetaList,
                Surfaces,
                encl,
                setting,
                options,
                tolerances,
                numeric_format,
            )
            voidList.append(voids)

    voidList.append(set_graveyard_cell(Surfaces, UniverseBox, options, tolerances, numeric_format))

    return VF.update_void_list(init, voidList, NestedEnclosure, setting.sort_enclosure)


def get_void_def(
    MetaList,
    Surfaces,
    Enclosure,
    setting,
    options,
    tolerances,
    numeric_format,
    Lev0=False,
):

    maxsurf = setting.maxSurf
    maxbracket = setting.maxBracket
    minSize = setting.minVoidSize

    if "full" in setting.simplify.lower():
        simplifyVoid = "full"
    elif "void" in setting.simplify.lower():
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
        logger.info("Loop, Box to Split :{iloop}, {nvoid}")

        for iz, z in enumerate(tqdm(Initial, desc=f"Void Generation Loop: {iloop}")):
            nsurfaces, nbrackets = z.get_numbers()
            logger.info(f"{iloop} {iz + 1}/{nvoid} {nsurfaces} {nbrackets}")

            if nsurfaces > maxsurf and nbrackets > maxbracket:
                newspace = z.split(minSize)
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

                logger.info(f"build complementary {iloop} {iz}")

                cell, CellIn = z.get_void_complementary(Surfaces, options, tolerances, numeric_format, simplify=simplifyVoid)
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
        mVoid.set_definition(vcell[0], simplify=True)
        mVoid.set_material(Enclosure.Material, Enclosure.Rho, Enclosure.MatInfo)
        mVoid.set_dilution(Enclosure.Dilution)

        mVoid.__commentInfo__ = vcell[1]

        voidList.append(mVoid)

    return voidList


def set_graveyard_cell(Surfaces, UniverseBox, options, tolerances, numeric_format):
    Universe = VoidBox([], UniverseBox)

    externalBox = get_universe_complementary(Universe, Surfaces, options, tolerances, numeric_format)
    center = UniverseBox.Center
    radius = 0.51 * UniverseBox.DiagonalLength
    sphere = GeounedSurface(("Sphere", (center, radius)), UniverseBox)
    id, exist = Surfaces.add_sphere(sphere, tolerances)

    sphdef = BoolSequence(str(-id))
    sphdef.operator = "AND"
    sphdef.append(externalBox)

    notsph = BoolSequence(str(id))

    mVoidSphIn = GeounedSolid(0)
    mVoidSphIn.Void = True
    mVoidSphIn.CellType = "void"
    mVoidSphIn.set_definition(sphdef)
    mVoidSphIn.set_material(0, 0, "Graveyard_in")
    mVoidSphIn.__commentInfo__ = None

    mVoidSphOut = GeounedSolid(1)
    mVoidSphOut.Void = True
    mVoidSphOut.CellType = "void"
    mVoidSphOut.set_definition(notsph)
    mVoidSphOut.set_material(0, 0, "Graveyard")
    mVoidSphOut.__commentInfo__ = None

    return (mVoidSphIn, mVoidSphOut)


# TODO check this is being used
def get_universe_complementary(Universe, Surfaces, options, tolerances, numeric_format):
    Def = BoolSequence(operator="OR")
    for p in Universe.get_bound_planes():
        id, exist = Surfaces.add_plane(p, options, tolerances, numeric_format, False)
        if not exist:
            Def.elements.append(-id)
        else:
            s = Surfaces.get_surface(id)
            if is_opposite(p.Surf.Axis, s.Surf.Axis):
                Def.elements.append(id)
            else:
                Def.elements.append(-id)
    return Def


def void_comment_line(CellInfo):
    boxDef, cellIn = CellInfo
    cells = ", ".join(map(str, cellIn))
    box = ", ".join(f"{num:.3f}" for num in boxDef)
    line = f"Automatic Generated Void Cell. Enclosure({box})\n"
    line += f"Enclosed cells : ({cells})\n"
    return line
