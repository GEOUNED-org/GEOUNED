#
# Module to load a STEP file
#
import os
import re
import FreeCAD, Part
from FreeCAD import Import
from ..Utils import Functions as UF
from . import LoadFunctions as LF

# Paco mod
def extractMaterials(filename):
    rhoreal = []
    namestr = []
    mdict = {} #_ Material dictionary
    with open(filename,"rt") as file:
      for line in file:
        vals = line.split()
        if vals[0].startswith("#"): continue
        matlabel = int(vals[0])
        rhoreal  = -float(vals[1])
        matname  = " ".join(vals[2:])
        mdict[matlabel] = (rhoreal,matname)
    return mdict

def LoadCAD(filename,matfilename,defaultMat=[],compSolids=True):

    # Set document solid tree options when opening CAD differing from version 0.18
    if int(FreeCAD.Version()[1]) > 18 : LF.set_docOptions()

    CAD_simplificado_doc=FreeCAD.newDocument("CAD_simplificado")
    Import.insert(filename,"CAD_simplificado")

    if matfilename != '' :
       if os.path.exists(matfilename): 
          mdict = extractMaterials(matfilename)
       else:
          print('Material definition file {} does not exist.'.format(matfilename))
          mdict = {} 
    else:
       mdict = {}
 
    s=Part.Shape()
    s.read(filename)
    Solids=s.Solids
    MetaList = []
    for i,s in enumerate(Solids) :
      MetaList.append(UF.GEOUNED_Solid(i+1,s))

    isolid = 0
    missingMat = set()

    docObjects = CAD_simplificado_doc.Objects

    for elem in docObjects:
        if (elem.TypeId  == 'Part::Feature'):
           comment = LF.getCommentTree(elem)
           if not elem.Shape.Solids:
              print('Warning: Element {:} has no associated solid'.format(comment+'/'+elem.Label))
              continue
           else: 
              tempremat = None
              tempredil = None
            
              # MIO: lightly modification of label if required
              label = LF.GetLabel(elem.Label)
              comment=comment+'/'+label
              if elem.InList:
                 # MIO: lightly modification of label if required
                 label_inList = LF.GetLabel(elem.InList[0].Label)
                 enclLabel  = re.search('enclosure(?P<encl>[0-9]+)_(?P<parent>[0-9]+)_', label_inList)
                 if not enclLabel :
                     enclLabel  = re.search('enclosure(?P<encl>[0-9]+)_(?P<parent>[0-9]+)_', label)
                     
                 envelLabel = re.search('envelope(?P<env>[0-9]+)_(?P<parent>[0-9]+)_', label_inList)
                 if not envelLabel :
                     envelLabel = re.search('envelope(?P<env>[0-9]+)_(?P<parent>[0-9]+)_', label)

                 #tempremat = re.search("(m(?P<mat>\d+)_)",elem.Label)
                 #if not tempremat :
                 #    tempremat = re.search("(m(?P<mat>\d+)_)",elem.InList[0].Label)

                 #tempredil = re.search("(_d(?P<dil>\d*\.\d*)_)",elem.Label)
                 #if not tempredil :
                 #    tempredil = re.search("(_d(?P<dil>\d*\.\d*)_)",elem.InList[0].Label)
                     
                 # Paco modifications
                 # Search for material definition in tree
                 xelem = [elem]
                 while xelem and not tempremat:
                    # MIO: Modification of label if required
                    temp_label = LF.GetLabel(xelem[0].Label)
                    tempremat = re.search("_m(?P<mat>\d+)_","_"+temp_label )
                    xelem = xelem[0].InList

                 # Search for dilution definition in tree
                 xelem = [elem]
                 while xelem and not tempredil:
                    # MIO: Modification of label if required
                    temp_label = LF.GetLabel(xelem[0].Label)
                    tempredil = re.search("_d(?P<dil>\d*\.\d*)_",temp_label)
                    xelem = xelem[0].InList
                 # Paco end
              else:
                 enclLabel   = None
                 envelLabel  = None

            # compSolid Diferent solid of the same cell are stored in the same metaObject (compSolid)
            # enclosures and envelopes are always stored as compound
              if  compSolids or enclLabel or envelLabel :
                 
                 init = isolid
                 end  = isolid+len(elem.Shape.Solids)
                 LF.fuseMetaObj(MetaList,init,end)
                 nsolids = 1
              else:
                 nsolids = len(elem.Shape.Solids)
              
              for i in range(nsolids):
                  MetaList[isolid].setComments('{}{}'.format(comment,(i+1)))
                  MetaList[isolid].setCADSolid()
                
                  if tempremat:
                      matlabel = int(tempremat.group('mat'))
                      if matlabel in mdict.keys() :
                         MetaList[isolid].setMaterial(matlabel,mdict[matlabel][0],mdict[matlabel][1])
                      else:
                         if matlabel == 0 :
                            MetaList[isolid].setMaterial(matlabel,0,0)
                         else:
                            MetaList[isolid].setMaterial(matlabel,-100,'Missing material density information')
                            missingMat.add(matlabel)          
                  else:
                     #print('Warning : No material label associated to solid {}.\nDefault material used instead.'.format(comment)) 
                     if defaultMat:
                         MetaList[isolid].setMaterial(*defaultMat)
                  if tempredil:
                      MetaList[isolid].setDilution(float(tempredil.group('dil'))) 
 
                  if enclLabel is not None:
                      MetaList[isolid].EnclosureID = int(enclLabel.group('encl'))
                      MetaList[isolid].ParentEnclosureID = int(enclLabel.group('parent'))
                      MetaList[isolid].IsEnclosure = True
                      MetaList[isolid].CellType    = 'void'
                       
                  if envelLabel is not None:
                      MetaList[isolid].EnclosureID = int(envelLabel.group('env'))
                      MetaList[isolid].ParentEnclosureID = int(envelLabel.group('parent'))
                      MetaList[isolid].IsEnclosure = True
                      MetaList[isolid].CellType    = 'envelope'  
                  isolid += 1     

    LF.joinEnvelopes(MetaList)
    if missingMat :
        print('Warning!! At least one material in the CAD model is not present in the material file')
        print('List of not present materials:', missingMat)

    EnclosureList = LF.setEnclosureSolidList(MetaList)
    if EnclosureList:
        LF.checkEnclosure(CAD_simplificado_doc,EnclosureList)
        # LF.RemoveEnclosure(MetaList)
        return MetaList,EnclosureList
    else:
        return MetaList,[]










