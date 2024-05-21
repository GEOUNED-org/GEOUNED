def assignEnclosure(MetaList, NestedEnclosure):
    """Assign to all solids the Enclosures ID it belongs to"""

    for m in MetaList:
        if m.IsEnclosure:
            continue
        dep = set()
        outside = True
        for lev, Level in enumerate(reversed(NestedEnclosure)):
            for encl in Level:
                inter = m.check_intersection(encl.CADSolid)
                # CAD solid intersect level i enclosure -> end of nested loops
                # Current loop continue because solid can intersect various enclosure
                if inter == 0:
                    dep.add(encl.EnclosureID)
                    m.ParentEnclosureID = encl.EnclosureID
                    if lev == len(NestedEnclosure) - 1:
                        dep.add(-1)
                        m.ParentEnclosureID = -1

                # CAD solid inside level i enclosure -> end of current loop go up i-1 level
                elif inter == -1:
                    dep.add(encl.EnclosureID)
                    m.ParentEnclosureID = encl.EnclosureID
                    outside = False
                    break

            # stop loop if CAD solid outside all level i enclosure or intersecting one
            if not outside:
                m.enclosure_list = dep
                break

        if not dep:
            dep.add(-1)
            m.ParentEnclosureID = -1
        if outside:
            m.enclosure_list = dep
    return


def select_solids(MetaList, LowLevelEnclosure, Enclosure):
    """This function selects the CAD solids and nested enclosures
    that can intersect with the current nested enclosure.
    For intermediate Levels selected solids are those intersecting i and i+1 level enclosure,
    or fully inside level i enclosure and fully outside i+1 enclosure"""

    newMetaList = LowLevelEnclosure[:]

    if "BoundBox" in str(Enclosure):
        enclID = -1
    else:
        enclID = Enclosure.EnclosureID

    for m in MetaList:
        if m.IsEnclosure:
            continue
        if enclID in m.enclosure_list:
            newMetaList.append(m)

    return newMetaList


def update_void_list(offset, voidList, NestedEnclosure, sort_enclosure):

    newVoidList = []
    if NestedEnclosure:
        update_comment = True
    else:
        update_comment = False

    icount = offset + 1
    voids = voidList[0]
    for m in voids:
        m.__id__ = icount
        if update_comment and not sort_enclosure:
            m.Comments = m.Comments + "\nLevel 0 void enclosure"
        icount += 1
        m.ParentEnclosureID = -1
        m.EnclosureID = 0
        newVoidList.append(m)

    ivoid = 1
    for Level in NestedEnclosure:
        for encl in Level:
            if encl.CellType == "envelope":
                continue
            for m in voidList[ivoid]:
                m.__id__ = icount
                m.ParentEnclosureID = encl.ParentEnclosureID
                m.EnclosureID = encl.EnclosureID
                if sort_enclosure:
                    m.Comments = m.Comments + f"\n{encl.Comments}"
                else:
                    m.Comments = m.Comments + f"\nVoid Enclosure #{encl.EnclosureID}"
                icount += 1
                newVoidList.append(m)
            ivoid += 1

    # graveyard Cell

    Graveyard = voidList[-1]
    for m in Graveyard:
        m.__id__ = icount
        m.ParentEnclosureID = -1
        m.EnclosureID = 0
        newVoidList.append(m)
        icount += 1

    return newVoidList
