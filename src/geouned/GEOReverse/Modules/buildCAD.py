import BOPTools.SplitAPI

from .buildSolidCell import FuseSolid
from .Utils.booleanFunction import BoolSequence


def buildCAD(UnivCell, data, config):

    UniverseCut = True
    if "Ustart" not in config.keys():
        config["Ustart"] = 0
    if "levelMax" not in config.keys():
        config["levelMax"] = "all"
    UnivCell.name = 0
    UnivCell.Fill = config["Ustart"]

    # read all surfaces definition
    if config["format"] == "mcnp":
        factor = 10
    else:
        factor = 1

    modelSurfaces = data.GetSurfaces(scale=factor)  # scale change cm in mcnp to mm in CAD Obj

    # read Cells and group into universes
    print(config)
    levels, UniverseCells, modelSurfaces = data.GetFilteredCells(modelSurfaces, config)

    # assign to each cell the surfaces belonging to the cell
    AssignSurfaceToCell(UniverseCells, modelSurfaces)

    #    print(UniverseCells[0][120].definition.str)
    #    print(UniverseCells[0][120].surfaces)
    #    CT=build_c_table_from_solids(UnivCell.shape,UniverseCells[0][70],option='full')
    #    print(CT)
    #    simply = BoolSequence(UniverseCells[0][70].definition.str)
    #    print('antesSimply',simply)
    #    simply.simplify(CT)
    #    print('despues',simply)
    # exit()

    # dictionnary of the cells filled with a given Universe U
    # universeContainers = get_universe_containers(levels,UniverseCells)

    UnivCell.level = None
    levelMax = config["levelMax"]
    Ustart = config["Ustart"]
    if levelMax == "all":
        levelMax = len(levels)

    for lev, Univ in levels.items():
        if Ustart in Univ:
            UnivCell.level = lev - 1
            break
    startInfo = (Ustart, levelMax)

    return BuildUniverse(startInfo, UnivCell, UniverseCells, universeCut=UniverseCut)


def interferencia(container, cell, mode="slice"):

    if mode == "common":
        return cell.shape.common(container.shape)

    Base = cell.shape
    Tool = (container.shape,)

    solids = BOPTools.SplitAPI.slice(Base, Tool, "Split", tolerance=0).Solids
    cellParts = []
    for s in solids:
        if container.shape.isInside(s.CenterOfMass, 0.0, False):
            cellParts.append(s)

    if not cellParts:
        return cell.shape
    else:
        return FuseSolid(cellParts)


def AssignSurfaceToCell(UniverseCells, modelSurfaces):

    for Uid, uniCells in UniverseCells.items():
        for cId, c in uniCells.items():
            c.setSurfaces(modelSurfaces)


def get_universe_containers(levels, Universes):
    Ucontainer = {}
    for lev in range(1, len(levels)):
        for U, name in levels[lev]:
            UFILL = Universes[U][name].FILL
            if UFILL in Ucontainer.keys():
                Ucontainer[UFILL].append((U, name, lev))
            else:
                Ucontainer[UFILL] = [(U, name, lev)]
    return Ucontainer


def BuildUniverse(startInfo, ContainerCell, AllUniverses, universeCut=True, duplicate=False):
    CADUniverse = []

    Ustart, levelMax = startInfo
    Universe = AllUniverses[Ustart]

    print(f"Build Universe {ContainerCell.FILL} in container cell {ContainerCell.name}")
    fails = []
    for NTcell in Universe.values():
        if duplicate:
            if NTcell.shape:
                buildShape = False
                if ContainerCell.CurrentTR:
                    cell = NTcell.copy()
                    cell.transformSolid(ContainerCell.CurrentTR)
                else:
                    cell = NTcell
            else:
                CTRF = None
                buildShape = True
                if ContainerCell.CurrentTR:
                    CC = ContainerCell.copy()
                    CC.transformSolid(CC.CurrentTR, reverse=True)
                else:
                    CC = ContainerCell
        else:
            CTRF = ContainerCell.CurrentTR
            CC = ContainerCell
            buildShape = True
            NTcell = NTcell.copy()

        if buildShape:
            print(f"Level :{CC.level + 1}  build Cell {NTcell.name} ")
            if type(NTcell.definition) is not BoolSequence:
                NTcell.definition = BoolSequence(NTcell.definition.str)

            bBox = ContainerCell.shape.BoundBox

            debug = False
            if debug:
                NTcell.buildShape(bBox, surfTR=CTRF, simplify=False)
            else:
                try:
                    NTcell.buildShape(bBox, surfTR=CTRF, simplify=False)
                except:
                    print(f"fail converting cell {NTcell.name}")
                    fails.append(NTcell.name)

            if NTcell.shape is None:
                continue

            if duplicate:
                if ContainerCell.CurrentTR:
                    cell = NTcell.copy()
                    cell.transformSolid(ContainerCell.CurrentTR)
                else:
                    cell = NTcell
            else:
                cell = NTcell

        if universeCut and ContainerCell.shape:
            cell.shape = interferencia(ContainerCell, cell)

        if not cell.FILL or ContainerCell.level + 1 == levelMax:
            CADUniverse.append(cell)
        else:
            if ContainerCell.CurrentTR:
                cell.CurrentTR = ContainerCell.CurrentTR.multiply(cell.TRFL)
            cell.level = ContainerCell.level + 1
            univ, ff = BuildUniverse((cell.FILL, levelMax), cell, AllUniverses, universeCut=universeCut)
            CADUniverse.append(univ)
            fails.extend(ff)

    return ((ContainerCell.name, Ustart), CADUniverse), fails


def makeTree(CADdoc, CADCells):

    label, universeCADCells = CADCells
    groupObj = CADdoc.addObject("App::Part", "Materials")

    groupObj.Label = f"Universe_{label[1]}_Container_{label[0]}"

    CADObj = {}
    for i, c in enumerate(universeCADCells):
        if type(c) is tuple:
            groupObj.addObject(makeTree(CADdoc, c))
        else:
            featObj = CADdoc.addObject("Part::FeaturePython", f"solid{i}")
            featObj.Label = f"Cell_{c.name}_{c.MAT}"
            featObj.Shape = c.shape
            if c.MAT not in CADObj.keys():
                CADObj[c.MAT] = [featObj]
            else:
                CADObj[c.MAT].append(featObj)

    for mat, matGroup in CADObj.items():
        groupMatObj = CADdoc.addObject("App::Part", "Materials")
        groupMatObj.Label = f"Material_{mat}"
        groupMatObj.addObjects(matGroup)
        groupObj.addObject(groupMatObj)

    return groupObj
