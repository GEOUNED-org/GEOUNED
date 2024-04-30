import re

import FreeCAD

from ..Utils.Functions import GeounedSolid
from ..Utils.Options.Classes import Options as opt


def GetLabel(label):
    """Deleting the last word of the string if this is a number
    Only if the option delLastNumber is True"""
    if not opt.delLastNumber:
        return label
    wrd = label.split()
    try:
        new_label = " ".join(wrd[:-1])
        return new_label
    except:
        return label


def getCommentTree(obj):
    recursive_list = []
    c_obj = obj
    while c_obj.InList:
        label = GetLabel(c_obj.InList[0].Label)
        recursive_list.append(label)
        c_obj = c_obj.InList[0]

    comment = ""
    for label in reversed(recursive_list):
        comment += "/" + label
    return comment


def joinEnvelopes(meta_list):
    join_list = []
    old_id = -1
    new_list = []
    for i, m in enumerate(meta_list):
        if m.CellType != "envelope":
            continue
        if m.EnclosureID != old_id:
            join_list.append(new_list)
            new_list = [i]
            old_id = m.EnclosureID
        else:
            new_list.append(i)

    join_list.append(new_list)
    del join_list[0]

    for id_list in join_list:
        fuseMetaObj(meta_list, id_list[0], id_list[-1] + 1)


def fuseMetaObj(meta_list, init, end):
    if end - init == 1:
        meta_list[init].__id__ = init + 1
        return
    solids = []
    for m in meta_list[init:end]:
        solids.extend(m.Solids)

    new_meta = GeounedSolid(init + 1, solids)
    new_meta.EnclosureID = meta_list[init].EnclosureID
    new_meta.ParentEnclosureID = meta_list[init].ParentEnclosureID
    new_meta.IsEnclosure = meta_list[init].IsEnclosure
    new_meta.CellType = meta_list[init].CellType

    del meta_list[init:end]
    meta_list.insert(init, new_meta)


# Paco mod
def checkCADFileMaterial(material, mdict):
    """Function to check that if a component in the CAD file has a material, it exist in the material file"""
    templist = [elem for elem in set(material) if elem not in mdict]
    if templist:
        raise ValueError(
            "At least one material in the CAD model is not present in the material file\n"
            f"List of not present materials:{templist}\n"
            "Code STOPS\n"
        )
    return


# Paco mod
def setEnclosureSolidList(meta_list):
    """Function to define DataList list."""
    return [elem for elem in meta_list if elem.IsEnclosure]


def RemoveEnclosure(meta_list):
    """This function removes Enclosure solid and actualises __id__ attributes to take into account the removal from meta_list list of nested enclosures originally loaded from CAD model."""

    # update id of solids
    len_meta_list = len(meta_list)
    i_count = 0
    for i, m in enumerate(meta_list):
        if m.IsEnclosure:
            i_count += 1
        else:
            meta_list[i].__id__ -= i_count

    # remove Enclosure solids from metaList
    for i in range(len_meta_list - 1, -1, -1):
        if meta_list[i].IsEnclosure:
            meta_list.pop(i)


def setEnclosureLevels(enclosure_list):
    nested_list = [[0], []]
    nested_enclosure = [[]]
    parent_level = 0
    temp_list = enclosure_list[:]

    while len(temp_list) > 0:
        remove = []
        for i, encl in enumerate(temp_list):
            if encl.ParentEnclosureID in nested_list[parent_level]:
                nested_list[parent_level + 1].append(encl.EnclosureID)
                nested_enclosure[parent_level].append(encl)
                remove.append(i)

        nested_list.append([])
        nested_enclosure.append([])
        parent_level += 1
        remove.reverse()
        for i in remove:
            del temp_list[i]

    nested_enclosure_list = []
    for Level in nested_enclosure:
        temp = []
        for i, encl in enumerate(Level):
            temp.append((encl.ParentEnclosureID, i))
        temp.sort()

        grouped = []
        for t in temp:
            grouped.append(Level[t[1]])
        nested_enclosure_list.append(grouped)

    for i, Level in enumerate(nested_enclosure_list[0:-1]):
        for encl in Level:
            child_list = []
            for child in nested_enclosure_list[i + 1]:
                if child.ParentEnclosureID == encl.EnclosureID:
                    child_list.append(child)
            encl.SonEnclosures = child_list[:]

    return nested_enclosure_list


def checkEnclosure(freecad_doc, enclosure_list):

    stop = False
    # check all enclosure labels have an associated solid
    temp_list = []
    for elem in freecad_doc.Objects:
        if elem.TypeId == "Part::Feature":
            if elem.Shape.Solids:
                if elem.InList:
                    templabel = re.search(
                        "enclosure(?P<encl>[0-9]+)_(?P<parent>[0-9]+)_",
                        elem.InList[0].Label,
                    )
                    if templabel is not None:
                        if (
                            elem.TypeId == "Part::Feature"
                            and len(elem.Shape.Solids) == 0
                        ):
                            temp_list.append(elem)

    if temp_list:
        stop = True
        print(
            "One or more nested enclosure labels in CAD solid tree view/structure tree do not have any CAD solid.\n",
            "Each nested enclosure must have only one solid.\nCode STOPS.",
            "\nList of problematic nested enclosure labels:",
        )
        for elem in temp_list:
            print(elem.EnclosureID)

    # check enclosure Labels don't make loops

    # 1) check at leat one 0 is in Parent enclosure label
    # 2) check 0 is not in child enclosure label
    # 3) check child enclosure label is not repeated

    sid_list = []
    pid_set = set()
    repeated_id = set()

    for encl in enclosure_list:
        pid_set.add(encl.ParentEnclosureID)
        if encl.EnclosureID in sid_list:
            repeated_id.add(encl.EnclosureID)
        else:
            sid_list.append(encl.EnclosureID)

    if 0 in sid_list:
        stop = True
        print('"0" cannot be label on child Enclosure')
    if 0 not in pid_set:
        stop = True
        print(' "0" should parent label of most external enclosure(s)')
    if repeated_id:
        stop = True
        print("Child label cannot be repeated.\nRepeated labels :")
        for lab in repeated_id:
            print(lab)

    # this stop should not be move to the end of routine
    # if previous point not satisfied point 4) may lead to infinite loop
    if stop:
        raise ValueError("Exiting to avoid infinite loop")

    # 4) look for explicit loops
    encl_dict = dict()
    for encl in enclosure_list:
        if encl.ParentEnclosureID in encl_dict.keys():
            encl_dict[encl.ParentEnclosureID].append(encl.EnclosureID)
        else:
            encl_dict[encl.ParentEnclosureID] = [encl.EnclosureID]

    parent = [0]
    cont = True
    while cont:
        cont = False
        next_parent = []
        for p in parent:
            if p in encl_dict.keys():
                cont = True
                next_parent.extend(encl_dict[p])
                del encl_dict[p]
        parent = next_parent

    if encl_dict.keys():
        print("Following enclosure produce loop")
        for p in encl_dict.keys():
            for s in encl_dict[p]:
                print(f"{s}_{p}")
        raise ValueError("GEOUNED.LoadFunctions.checkEnclosure failed")

    # check enclosures solids are correctly nested and don't overlap each other
    nested_levels = setEnclosureLevels(enclosure_list)

    overlap = []
    encl_tree = [[0]]
    for level in nested_levels:
        encl_tree = updateTree(encl_tree, level)
        for encl in level:
            same_parent = dict()
            if encl.ParentEnclosureID in same_parent.keys():
                same_parent[encl.ParentEnclosureID].append(encl.EnclosureID)
            else:
                same_parent[encl.ParentEnclosureID] = [encl.EnclosureID]
        for encl in same_parent.values():
            overlap.extend(checkOverlap(encl))

    not_embedded = []
    for chain in encl_tree:
        up = chain[0]
        for low in chain[1:]:
            inter = up.checkIntersection(low.CADSolid)
            if inter != -2:
                not_embedded.append((low.EnclosureID, up.EnclosureID))

    if not_embedded:
        stop = True
        print(" Following enclosures are not fully embedded in Parent enclosure")
        for elemt in not_embedded:
            print("{}_{}").format(elemt[0], elemt[1])

    if overlap:
        stop = True
        print(" Following enclosures overlapping ")
        for elemt in overlap:
            print("{}_{}").format(elemt[0], elemt[1])

    if stop:
        raise ValueError("GEOUNED.LoadFunctions.checkEnclosure failed")


def checkOverlap(enclosures):
    overlap = []
    for i, enc1 in enumerate(enclosures):
        for enc2 in enclosures[i + 1 :]:
            inter = enc1.checkIntersection(enc2.CADSolid)
            if inter != 1:
                overlap.append((enc1.EnclosureID, enc2.EnclosureID))
    return overlap


def updateTree(Tree, level):
    new_tree = []
    for encl in level:
        for lst in Tree:
            if lst[-1] == encl.ParentEnclosureID:
                new_tree.append(lst + [encl.EnclosureID])
                continue
    return new_tree


def set_docOptions():
    # set import step document options for FreeCAD version >0.18 compatible with FreeCAD0.18 opening options
    p0 = FreeCAD.ParamGet("User parameter:BaseApp/Preferences/Mod/Import")
    p0.SetBool("UseLinkGroup", False)
    p0.SetBool("ReduceObjects", False)
    p0.SetBool("ExpandCompound", False)

    p1 = FreeCAD.ParamGet("User parameter:BaseApp/Preferences/Mod/Import/hSTEP")
    p1.SetBool("UseLinkGroup", False)
    p1.SetBool("ReadShapeCompoundMode", True)


# functions not used in the code
def checkIndex(docList, index, returnObject=False):
    subList = docList.ElementList
    for i in index[:-1]:
        if len(subList) > i:
            subList = subList[i]
            if subList.TypeId != "App::LinkGroup":
                return None
            subList = subList.ElementList
        else:
            return None

    if len(subList) > index[-1]:
        if subList[index[-1]].TypeId == "App::LinkGroup":
            return "Link"
        else:
            if returnObject:
                return subList[index[-1]]
            else:
                return "Feature"
    else:
        return None


def nextIndex(docList, lastIndex=None):
    if lastIndex is None:
        lastIndex = [0]
        while checkIndex(docList, lastIndex) != "Feature":
            lastIndex.append(0)
        else:
            return lastIndex

    tmp = lastIndex[:]
    move = True
    while True:
        if move:
            tmp[-1] = tmp[-1] + 1
        obj = checkIndex(docList, tmp)
        if obj == "Feature":
            return tmp
        elif obj == "Link":
            tmp.append(0)
            move = False
        else:
            tmp = tmp[0:-1]
            move = True
            if len(tmp) == 0:
                return None


def solidGenerator(doclist):
    last = None
    while True:
        last = nextIndex(doclist, last)
        if last is None:
            break
        yield checkIndex(doclist, last, True)
