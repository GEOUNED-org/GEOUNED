import Part

from .options import Options
from .splitFunction import SplitBase, SplitSolid, joinBase
from .Utils.booleanFunction import BoolSequence


def getPart(slist):
    sol = []
    for s in slist:
        if type(s) is list:
            sol.extend(getPart(s))
        else:
            sol.append(s)
    return sol


def BuildSolid(cell, boundBox, mode="oneByOne", simplify=False):

    cutCell = cell.makeBox(boundBox)
    # cell.definition = BoolSequence(cell.definition.str)
    cell.cleanUndefined()

    celParts = BuildDepth(cell, SplitBase(cutCell), mode, True, simplify)

    celParts = getPart(celParts)
    # print('celparts',len(celParts))
    shapeParts = []
    for i, s in enumerate(celParts):
        shapeParts.append(s.base)
        # s.base.exportStep('solid{}.stp'.format(i))

    #   tt = FuseSolid(shapeParts)
    #   tt = tt.removeSplitter()
    # tt=Part.makeCompound(shapeParts)
    # tt.exportStep('cell{}.stp'.format(cell.name))
    return shapeParts
    # return FuseSolid(shapeParts)


def BuildDepth(cell, cutShape, mode, baseBox, simplify=False, loop=0):

    loop += 1
    seq = cell.definition
    if seq.level == 0:
        if baseBox:
            cutShape, cut = BuildSolidParts(cell, cutShape, mode)
        else:
            cutShape, cut = BuildSolidParts(cell, cutShape, "solids")
        return cutShape

    if type(cutShape) is not list:
        cutShape = [cutShape]
    newCutShape = []
    for i, CS in enumerate(cutShape):
        cbaseBox = baseBox
        # CS.base.exportStep('CS_{}_{}.stp'.format(i,str(cell.definition)))
        # CTable =build_c_table_from_solids(cell.makeBox(CS.base.BoundBox),cell.surfaces,option='full')
        # cell.definition.simplify(CTable)
        cell.definition.group_single()

        if type(cell.definition.elements) is not bool:
            if cell.definition.level == 0:
                tmp = BoolSequence(operator=cell.definition.operator)
                tmp.append(cell.definition)
                cell.definition = tmp

            if seq.operator == "AND":
                part = CS
                for e in cell.definition.elements:
                    part = BuildDepth(cell.getSubCell(e), part, mode, cbaseBox, simplify, loop=loop)
                    cbaseBox = False
                newCutShape.extend(part)
            else:
                cellParts = []
                for e in cell.definition.elements:
                    sub = cell.getSubCell(e)
                    part = BuildDepth(sub, CS, mode, baseBox, simplify, loop=loop)
                    cellParts.extend(part)

                JB = joinBase(cellParts)
                if JB.base is not None:
                    newCutShape.append(JB)

        elif cell.definition.elements:
            newCutShape.append(CS)

    cutShape = newCutShape
    return cutShape


def BuildSolidParts(cell, base, mode):

    # part if several base in input
    if type(base) is list or type(base) is tuple:
        fullPart = []
        cutPart = []
        for b in base:
            fullList, cutList = BuildSolidParts(cell, b, mode)
            fullPart.extend(fullList)
            cutPart.extend(cutList)
        return fullPart, cutPart

    boundBox = base.base.BoundBox
    surfaces = tuple(cell.surfaces.values())
    # print('\nbuild Parts :',mode)
    # print(cell.definition)
    # print(boundBox)

    if mode == "solids":

        # TODO consider making this buildShape call conditional
        # if cell.definition.operator == "OR" and False:
        #     Def = cell.definition
        #     cell.definition = cell.definition.get_complementary()
        #     cell.buildShape(boundBox, force=False, simplify=False)
        #     cell.definition = Def
        # else:
        #     cell.buildShape(boundBox, force=True, simplify=False, fuse=True)
        cell.buildShape(boundBox, force=True, simplify=False, fuse=True)

        # print('export')
        # base.base.exportStep('base.stp')
        # name=''.join(str(cell.definition))
        # cell.shape.exportStep('sol{}.stp'.format(name))
    else:
        cell.buildSurfaceShape(boundBox)

    if not surfaces:
        print("not cutting surfaces")
        return tuple(base.base), tuple()
    if mode == "solids":
        full, cut = SplitSolid(base, surfaces, cell, solidTool=True, tolerance=Options.splitTolerance)
    elif mode == "allSurfaces":
        full, cut = SplitSolid(base, surfaces, cell, tolerance=Options.splitTolerance)

    elif mode == "planeFirst":
        planes = []
        others = []
        for s in surfaces:
            # s.buildShape( boundBox)
            # s.shape.exportStep('Tool_{}_{}.stp'.format(s.type,s.id))
            if s.type == "plane":
                planes.append(s)
            else:
                others.append(s)

        if planes:

            full, cut = SplitSolid(base, planes, cell, tolerance=Options.splitTolerance)
            # for i,s in enumerate(full):
            #    s.exportStep('fullplane_{}.stp'.format(i))
            # for i,s in enumerate(cut):
            #    s.base.exportStep('cutplane_{}.stp'.format(i))
            # print('planes',full)
            # print('planes',cut)
        else:
            full = []
            cut = base

        if others:
            newf, cut = SplitSolid(cut, others, cell, tolerance=Options.splitTolerance)
            # print('others',newf)
            # print('others',cut)
        else:
            newf = []
        full.extend(newf)

    elif mode == "otherFirst":
        planes = []
        others = []
        for s in surfaces:
            # s.buildShape( boundBox)
            if s.type == "plane":
                planes.append(s)
            else:
                others.append(s)

        if others:
            full, cut = SplitSolid(base, others, cell, tolerance=Options.splitTolerance)
            # print('others',full)
            # print('others',cut)
        else:
            full = []
            cut = base

        if planes:
            newf, cut = SplitSolid(cut, planes, cell, tolerance=Options.splitTolerance)
            # print('planes',newf)
            # print('planes',cut)

        else:
            newf = []

        full.extend(newf)

    elif mode == "oneByOne":
        planes = []
        others = []
        for s in surfaces:
            if s.type == "plane":
                planes.append(s)
            else:
                others.append(s)

        if planes:
            full, cut = SplitSolid(base, planes, cell, tolerance=Options.splitTolerance)
        else:
            full = []
            cut = base

        # cut[0].base.exportStep('cutPlane.stp')
        for surf in others:
            newf, cut = SplitSolid(cut, (surf,), cell, tolerance=Options.splitTolerance)
            full.extend(newf)

    elif mode == "otherOneByOne":
        planes = []
        others = []
        for s in surfaces:
            if s.type == "plane":
                planes.append(s)
            else:
                others.append(s)

        cut = base
        full = []
        for surf in others:
            newf, cut = SplitSolid(cut, (surf,), cell, tolerance=Options.splitTolerance)
            full.extend(newf)

        for surf in planes:
            newf, cut = SplitSolid(cut, (surf,), cell, tolerance=Options.splitTolerance)
            full.extend(newf)

    return full, cut


def FuseSolid(parts):
    if (len(parts)) <= 1:
        if parts:
            solid = parts[0]
        else:
            return None
    else:
        try:
            fused = parts[0].fuse(parts[1:])
        except:
            fused = None

        if fused is not None:
            try:
                refinedfused = fused.removeSplitter()
            except:
                refinedfused = fused

            if refinedfused.isValid():
                solid = refinedfused
            else:
                if fused.isValid():
                    solid = fused
                else:
                    solid = Part.makeCompound(parts)
        else:
            solid = Part.makeCompound(parts)

    if solid.Volume < 0:
        solid.reverse()
    return solid


def noOR(Seq):
    if len(Seq.elements) == 1:
        # Seq.operator = 'AND'
        return Seq
    newOR = BoolSequence(operator="OR")
    neg = []
    for e in Seq.elements:
        AND = BoolSequence(operator="AND")
        AND.append(*neg, e)
        newOR.append(AND)
        neg.append(-e)
    return newOR
