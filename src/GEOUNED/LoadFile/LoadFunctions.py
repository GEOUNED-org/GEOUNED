import re
import FreeCAD
from ..Utils.Functions import GEOUNED_Solid
from ..Utils.Options.Classes import Options as opt


def GetLabel(label):
    """ Deleting the last word of the string if this is a number
        Only if the option delLastNumber is True"""
    if not opt.delLastNumber: return label
    wrd = label.split()
    try:
        val=float(wrd[-1])
        new_label = ' '.join(wrd[:-1])
        return new_label
    except: return label

def getCommentTree(obj):
    recursiveList = []
    cObj = obj
    while cObj.InList :
       label = GetLabel(cObj.InList[0].Label)
       recursiveList.append( label )
       cObj = cObj.InList[0]

    comment = ''
    for label in reversed(recursiveList):   
       comment += '/' + label
    return comment 

def joinEnvelopes(MetaList):
   joinList=[]
   oldID = -1
   newlist=[]
   for i,m in enumerate(MetaList):
      if m.CellType != 'envelope' : continue
      if m.EnclosureID != oldID:
         joinList.append(newlist)
         newlist = [i]
         oldID = m.EnclosureID
      else:   
         newlist.append(i)

   joinList.append(newlist)
   del joinList[0]

   for IDlist in joinList:
      fuseMetaObj(MetaList,IDlist[0],IDlist[-1]+1)
      
def fuseMetaObj(MetaList,init,end):
   if end-init == 1 :
      MetaList[init].__id__ = init+1
      return
   solids=[]
   for m in MetaList[init:end]:
      solids.extend(m.Solids)

   newMeta = GEOUNED_Solid(init+1,solids)
   newMeta.EnclosureID = MetaList[init].EnclosureID
   newMeta.ParentEnclosureID = MetaList[init].ParentEnclosureID
   newMeta.IsEnclosure = MetaList[init].IsEnclosure
   newMeta.CellType    = MetaList[init].CellType

   del MetaList[init:end]
   MetaList.insert(init,newMeta)
    

# Paco mod
def checkCADFileMaterial(Material, mdict):
    """Function to check that if a component in the CAD file has a material, it exist in the material file"""
    templist = [ elem for elem in set(Material) if elem not in mdict ]
    if templist:
        raise ValueError(
         'At least one material in the CAD model is not present in the material file\n'
         f'List of not present materials:{templist}\n' 
         'Code STOPS\n'
        )
    return

# Paco mod
def setEnclosureSolidList(MetaList):
    """Function to define DataList list."""
    return [ elem for elem in MetaList if elem.IsEnclosure ]

def RemoveEnclosure(MetaList):
    """This function removes Enclosure solid and actualises __id__ attributes to take into account the removal from MetaList list of nested enclosures originally loaded from CAD model."""
    
    # update id of solids
    lenMetaList = len(MetaList)
    icount = 0
    for i, m in enumerate(MetaList):
        if m.IsEnclosure:
            icount += 1
        else:
            MetaList[i].__id__ -= icount
            
    # remove Enclosure solids from metaList
    for i in range(lenMetaList-1,-1,-1):
        if MetaList[i].IsEnclosure :
            MetaList.pop(i)
                    

def setEnclosureLevels(EnclosureList):
    nestedList      = [[0],[]]
    nestedEnclosure = [[]]
    parentLevel = 0
    tempList = EnclosureList[:]
    
    while len(tempList) > 0 :
       remove = [] 
       for i,encl in enumerate(tempList):
          if encl.ParentEnclosureID in nestedList[parentLevel] :
              nestedList[parentLevel+1].append(encl.EnclosureID)
              nestedEnclosure[parentLevel].append(encl)
              remove.append(i)
             
       nestedList.append([])
       nestedEnclosure.append([])
       parentLevel += 1
       remove.reverse()
       for i in remove:
          del tempList[i]

    nestedEnclosureList = []     
    for Level in nestedEnclosure :
       temp=[]
       for i,encl in enumerate(Level):
          temp.append((encl.ParentEnclosureID,i))
       temp.sort()

       grouped = []
       for t in temp:
          grouped.append(Level[t[1]])
       nestedEnclosureList.append(grouped)

    for i,Level in enumerate(nestedEnclosureList[0:-1]) :
        for encl in Level:
           childList = []
           for child in nestedEnclosureList[i+1]:
              if child.ParentEnclosureID == encl.EnclosureID :
                  childList.append(child)
           encl.SonEnclosures = childList[:]       
    
    return  nestedEnclosureList   


def checkEnclosure(FreeCAD_doc,EnclosureList):

    stop = False
    # check all enclosure labels have an associated solid
    TempList = []
    for elem in FreeCAD_doc.Objects:
       if (elem.TypeId=='Part::Feature'):
         if elem.Shape.Solids:
           if elem.InList:
              templabel = re.search('enclosure(?P<encl>[0-9]+)_(?P<parent>[0-9]+)_', elem.InList[0].Label)
              if templabel is not None:
                 if elem.TypeId == 'Part::Feature' and len(elem.Shape.Solids) == 0 :
                     TempList.append(elem)

    if TempList:
       stop = True
       print('One or more nested enclosure labels in CAD solid tree view/structure tree do not have any CAD solid.\n',\
              'Each nested enclosure must have only one solid.\nCode STOPS.',\
              '\nList of problematic nested enclosure labels:')
       for elem in TempList:
           print(elem.EnclosureID)
       
    # check enclosure Labels don't make loops

    # 1) check at leat one 0 is in Parent enclosure label
    # 2) check 0 is not in child enclosure label
    # 3) check child enclosure label is not repeated
    
    SIDList = []
    PIDSet  = set()
    repeatedID = set()

    for encl in EnclosureList:
       PIDSet.add(encl.ParentEnclosureID)
       if encl.EnclosureID in SIDList:
          repeatedID.add(encl.EnclosureID)
       else:       
          SIDList.append(encl.EnclosureID)

    if 0 in SIDList :
        stop = True
        print ('"0" cannot be label on child Enclosure')
    if 0 not in PIDSet :
        stop = True
        print (' "0" should parent label of most external enclosure(s)')
    if repeatedID:
        stop = True
        print('Child label cannot be repeated.\nRepeated labels :')
        for lab in repeatedID:
            print(lab)

    # this stop should not be move to the end of routine
    # if previous point not satisfied point 4) may lead to infinite loop
    if stop : raise ValueError('Exiting to avoid infinite loop')

    # 4) look for explicit loops
    enclDict = dict()
    for encl in EnclosureList:
        if encl.ParentEnclosureID in enclDict.keys():
            enclDict[encl.ParentEnclosureID].append(encl.EnclosureID)
        else:    
            enclDict[encl.ParentEnclosureID] = [encl.EnclosureID]

    parent = [0]
    cont = True
    while cont:
       cont = False 
       nextParent = [] 
       for p in parent:
          if p in enclDict.keys():
             cont = True
             nextParent.extend(enclDict[p])
             del enclDict[p]
       parent = nextParent   

    if enclDict.keys():
        print('Following enclosure produce loop')
        for p in enclDict.keys():
            for s in enclDict[p] :
                print ('{}_{}'.format(s,p))
        raise ValueError('GEOUNED.LoadFunctions.checkEnclosure failed')

    # check enclosures solids are correctly nested and don't overlap each other
    nestedLevels = setEnclosureLevels(EnclosureList)

    overlap = []
    enclTree = [[0]]                   
    for level in nestedLevels :
        enclTree = updateTree(enclTree,level)               
        for encl in level:
           sameParent= dict()
           if encl.ParentEnclosureID in sameParent.keys():
               sameParent[encl.ParentEnclosureID].append(encl.EnclosureID)
           else:            
               sameParent[encl.ParentEnclosureID] = [encl.EnclosureID]
        for encl in sameParent.values():
            overlap.extend(checkOverlap(encl))

    notEmbedded = []
    for chain in enclTree :
        up = chain[0]               
        for low in chain[1:]:
           inter = up.checkIntersection(low.CADSolid)
           if inter != -2 :
             notEmbedded.append((low.EnclosureID,up.EnclosureID))          

    if notEmbedded :
       stop = True
       print ( ' Following enclosures are not fully embedded in Parent enclosure')
       for elemt in notEmbedded:
           print('{}_{}').format(elemt[0],elemt[1])
                       
    if overlap :
       stop = True
       print ( ' Following enclosures overlapping ')
       for elemt in overlap:
           print('{}_{}').format(elemt[0],elemt[1])
 
    if stop : raise ValueError('GEOUNED.LoadFunctions.checkEnclosure failed')
                       
def checkOverlap(enclosures):
     overlap = []
     for i,enc1 in enumerate(enclosures):
        for enc2 in enclosures[i+1:]:
           inter = enc1.checkIntersection(enc2.CADSolid)
           if inter != 1:
              overlap.append((enc1.EnclosureID,enc2.EnclosureID))
     return overlap                  

def updateTree(Tree,level):
    newTree = []
    for encl in level :
       for lst in Tree:
           if lst[-1] == encl.ParentEnclosureID :
               newTree.append(lst + [encl.EnclosureID])
               continue
    return newTree                  

def set_docOptions():
    # set import step document options for FreeCAD version >0.18 compatible with FreeCAD0.18 opening options
    p0 = FreeCAD.ParamGet("User parameter:BaseApp/Preferences/Mod/Import")
    p0.SetBool("UseLinkGroup",False)
    p0.SetBool("ReduceObjects",False)
    p0.SetBool("ExpandCompound",False)
    
    p1 = FreeCAD.ParamGet("User parameter:BaseApp/Preferences/Mod/Import/hSTEP")
    p1.SetBool("UseLinkGroup",False)
    p1.SetBool("ReadShapeCompoundMode",True)
    

# functions not used in the code                       
def checkIndex(docList,index,returnObject=False):
    subList = docList.ElementList
    for i in index[:-1]:
        if len(subList) > i :
           subList = subList[i]
           if subList.TypeId != 'App::LinkGroup'   :  return None
           subList=subList.ElementList
        else:
           return None

    if len(subList) > index[-1] :
       if subList[index[-1]].TypeId == 'App::LinkGroup'  :  
          return 'Link'
       else:
          if returnObject :
             return subList[ index[-1]]
          else: 
             return 'Feature'       
    else:
       return None    
       
def nextIndex(docList,lastIndex=None):
    if lastIndex is None:
       lastIndex = [0]
       while checkIndex(docList,lastIndex) != 'Feature' :
           lastIndex.append(0)    
       else:
           return lastIndex         

    tmp = lastIndex[:]
    move = True
    while True:
       if move:
          tmp[-1] = tmp[-1]+1
       obj = checkIndex(docList,tmp)    
       if obj == 'Feature' :
          return tmp
       elif  obj == 'Link'  : 
          tmp.append(0)
          move = False
       else:
          tmp =tmp[0:-1]
          move = True
          if len(tmp) == 0 : return None


def solidGenerator(doclist) :
   last = None
   while True:
      last = nextIndex(doclist,last)
      if last is None : break
      yield checkIndex(doclist,last,True)
