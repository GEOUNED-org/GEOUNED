"""File with the VoidBox class"""
import FreeCAD
import Part
import GEOUNED.Decompose.Decom_one as Decom
import GEOUNED.Conversion.CellDefinition as Conv

from GEOUNED.Utils.Options.Classes import Options as opt
from GEOUNED.Utils.BasicFunctions_part1 import isOposite
from GEOUNED.Utils.Functions import GEOUNED_Solid, GEOUNED_Surface
from GEOUNED.Utils.booleanFunction import BoolSequence
from GEOUNED.Utils.BooleanSolids import buildCTableFromSolids, removeExtraSurfaces


class VoidBox (): 
  def __init__(self,MetaSolids,Enclosure):

    self.Objects =[]
    if 'BoundBox' in str(Enclosure) :
        self.BoundBox = Enclosure
        self.PieceEnclosure = None
    else:      
        self.BoundBox = Enclosure.BoundBox      
        self.PieceEnclosure = Enclosure

    for m in MetaSolids :
       if not m.BoundBox : continue
       if m.BoundBox.isValid() :
         if self.BoundBox.intersect(m.BoundBox) :
           Obj = self.__copyMeta__(m)
           self.__removeExtraComp__(Obj,self.BoundBox)
           self.Objects.append(Obj)    
    return
  
  def Split(self,minSize=200):
      
    dims = [self.BoundBox.XLength,self.BoundBox.YLength,self.BoundBox.ZLength]
    coord = ['X','Y','Z']
    for i in range(2,-1,-1):
      if 0.5*dims[i] < minSize :
        del(dims[i])
        del(coord[i])

    if len(dims) == 0 : return None   
      
    ip = dims.index(max(dims))
#    print ('dims : {} {}'.format(coord[ip],0.5*dims[ip]))
    if coord[ip] == 'X' : 
       pos   = self.BoundBox.XMin + 0.5*self.BoundBox.XLength
       X1Min = self.BoundBox.XMin
       X1Max = pos
       X2Min = pos
       X2Max = self.BoundBox.XMax
       Y1Min = self.BoundBox.YMin
       Y1Max = self.BoundBox.YMax
       Y2Min = Y1Min
       Y2Max = Y1Max
       Z1Min = self.BoundBox.ZMin
       Z1Max = self.BoundBox.ZMax
       Z2Min = Z1Min
       Z2Max = Z1Max
    elif coord[ip] == 'Y' :
       pos   = self.BoundBox.YMin + 0.5*self.BoundBox.YLength 
       X1Min = self.BoundBox.XMin
       X1Max = self.BoundBox.XMax
       X2Min = X1Min
       X2Max = X1Max
       Y1Min = self.BoundBox.YMin
       Y1Max = pos
       Y2Min = pos
       Y2Max = self.BoundBox.YMax
       Z1Min = self.BoundBox.ZMin
       Z1Max = self.BoundBox.ZMax
       Z2Min = Z1Min
       Z2Max = Z1Max  
    else  :
       pos   = self.BoundBox.ZMin + 0.5*self.BoundBox.ZLength 
       X1Min = self.BoundBox.XMin
       X1Max = self.BoundBox.XMax
       X2Min = X1Min
       X2Max = X1Max
       
       Y1Min = self.BoundBox.YMin
       Y1Max = self.BoundBox.YMax 
       Y2Min = Y1Min
       Y2Max = Y1Max
       
       Z1Min = self.BoundBox.ZMin
       Z1Max = pos
       Z2Min = pos
       Z2Max = self.BoundBox.ZMax

    VMin1 = FreeCAD.Vector(X1Min,Y1Min,Z1Min)
    VMax1 = FreeCAD.Vector(X1Max,Y1Max,Z1Max)
    VMin2 = FreeCAD.Vector(X2Min,Y2Min,Z2Min)
    VMax2 = FreeCAD.Vector(X2Max,Y2Max,Z2Max)
    box1 = FreeCAD.BoundBox(VMin1,VMax1)
    box2 = FreeCAD.BoundBox(VMin2,VMax2)
    
    if self.PieceEnclosure == None:
        Space1 = VoidBox(self.Objects,box1)
        Space2 = VoidBox(self.Objects,box2)
        VoidBoxTuple = (Space1,Space2)
    else:
        Space1 = self.PieceEnclosureSplit(box1)
        Space2 = self.PieceEnclosureSplit(box2)
        VoidBoxTuple = (Space1,Space2)
        if Space1 == None:
            VoidBoxTuple = (Space2,)
        if Space2 == None:
            VoidBoxTuple = (Space1,)
        if Space1 == None and Space2 == None:
            VoidBoxTuple = ()

    return VoidBoxTuple    


  def PieceEnclosureSplit(self,Box,Tolerance=1.0e-13):
    """This function creates a box-shaped solid with the new limits of given bounding box and
  it is intersected with the piece of nested enclosure to create the new void cell.
  If the limited region does not intersect with the piece, no void cell is created."""
    
    Cube = Part.makeBox(
        Box.XLength, Box.YLength, Box.ZLength,
        FreeCAD.Vector(Box.XMin, Box.YMin, Box.ZMin),
        FreeCAD.Vector(0,0,1))
    dist = Cube.distToShape(self.PieceEnclosure)[0]
    try:
        if abs(dist/Box.DiagonalLength) > Tolerance:
            return None
    except ZeroDivisionError:
        return None
    ShapeObject = Cube.common(self.PieceEnclosure)
    try:
        reldif = (Cube.Volume-ShapeObject.Volume)/Cube.Volume
    except ZeroDivisionError:
        return None
    if abs(reldif) <= Tolerance:
        return VoidBox(self.Objects,Box)
    elif ShapeObject.Solids :
        Solid = ShapeObject.Solids[0]
        return VoidBox(self.Objects,Solid)
    else:
        return None 

      
  def refine(self):
     Cube = Part.makeBox(self.BoundBox.XLength, self.BoundBox.YLength, self.BoundBox.ZLength,
                         FreeCAD.Vector(self.BoundBox.XMin, self.BoundBox.YMin, self.BoundBox.ZMin),
                         FreeCAD.Vector(0,0,1))
     refinedList = []
     newcom = []
       
     for m in self.Objects:
        self.__removeExtraComp__(m,Cube,mode='dist')
     return   

  def getVoidComplementary(self,Surfaces,simplify='No'):

    if self.PieceEnclosure is None:
        boxDef = BoolSequence(operator='AND')
        center = self.BoundBox.Center
        bBox = self.BoundBox
        for p in self.getBoundPlanes():
           id,exist = Surfaces.addPlane(p)
           if exist:
              s = Surfaces.getSurface(id)
              if isOposite(p.Surf.Axis,s.Surf.Axis):
                  id = -id
           if isOposite(p.Surf.Axis,p.Surf.Position-center):
               boxDef.elements.append(id)
           else:
               boxDef.elements.append(-id)
           enclosure = False
        d = opt.enlargeBox
        
    else:
        UniverseBox = self.PieceEnclosure.BoundBox
        TempPieceEnclosure = GEOUNED_Solid(None,self.PieceEnclosure)
        comsolid,err=Decom.SplitSolid(Part.makeCompound(TempPieceEnclosure.Solids),UniverseBox)
        Surfaces.extend(Decom.ExtractSurfaces(comsolid,'All',UniverseBox,MakeObj=True))
        TempPieceEnclosure.updateSolids(comsolid.Solids)
        Conv.cellDef(TempPieceEnclosure,Surfaces,UniverseBox)

        boxDef = TempPieceEnclosure.Definition
        bBox = self.PieceEnclosure.BoundBox
        enclosure = True
        d = max(opt.enlargeBox,2)

    Box = Part.makeBox(bBox.XLength+2*d, bBox.YLength+2*d, bBox.ZLength+2*d,
        FreeCAD.Vector(bBox.XMin-d, bBox.YMin-d, bBox.ZMin-d),
        FreeCAD.Vector(0,0,1))       

    voidSolidDef = BoolSequence(operator='OR')

    for m in self.Objects:
        voidSolidDef.append(m.Definition)

    voidSolidDef.joinOperators()

    if not voidSolidDef.elements : 
         return boxDef     #  Cell to get complementary are null => Void is only box definition
       
    # join all basic solids into one big meta Object
    # CAD solid representation is not needed because
    # here we are working with surfaces and void box


    complementary = BoolSequence(operator='AND')
    complementary.append(boxDef)

    if (simplify != 'no' ) :
       surfList = set(voidSolidDef.getSurfacesNumbers())

       if enclosure : 
           surfList.update(set(boxDef.getSurfacesNumbers()))
       else:
           for s in boxDef.elements:
              val = True if s >0 else False
              voidSolidDef.substitute(abs(s),val)
           res = voidSolidDef.clean()

       if enclosure or res is None :
          surfaceDict = {}
          for surf in Surfaces.values():
             for s in surf:
                if s.Index in surfList:
                   surfaceDict[s.Index] = s

          CTable = buildCTableFromSolids(Box,surfaceDict,option=simplify)
          voidSolidDef = removeExtraSurfaces(voidSolidDef,CTable) 
    
    if voidSolidDef.elements == True :
        return None
    elif voidSolidDef.elements == False or voidSolidDef.elements == []:
        return boxDef

    if voidSolidDef.level == 0:
       compSeq = voidSolidDef.getComplementary()
    else:
       if voidSolidDef.level == 1 and voidSolidDef.operator == 'AND':
           compSeq = BoolSequence(operator='OR')
       else:
           compSeq = BoolSequence(operator='AND')
       for comp in voidSolidDef.elements :
          if simplify == 'no':
             chk =  comp.check()
             #solid in cover full Void cell volume  => Void cell doesn't exist
             if chk == True :
                 if opt.verbose: print('warning void Cell should not exist')
                 return None

             #solid cell is not in void cell Void cell volume  => doesn't contribute to void definition
             elif chk == False :
                 continue

          pmoc = comp.getComplementary()
          compSeq.append(pmoc)
   
    if simplify == 'full':
       if enclosure :
          complementary.append(compSeq)
          complementary.simplify(CTable)
       else :  
          compSeq.simplify(CTable)
          complementary.append(compSeq)
    else:
       compSeq.simplify(None)
       complementary.simplify(None)
       complementary.append(compSeq)
 

    res = complementary.clean()
    complementary.levelUpdate()
    
    if type(res) is bool :
       return None
    else :
       return complementary

  def getBoxNumber(self):
      return len(self.Objects)

  def getNumbers(self):
      ns = 0
      nb = 0
            
      for m in self.Objects:
        ns += len(m.Surfaces)
        nb += len(m.Definition.elements)

      return ns,nb

  def getBoundPlanes(self):
     Xmid = 0.5*(self.BoundBox.XMin+self.BoundBox.XMax)
     Ymid = 0.5*(self.BoundBox.YMin+self.BoundBox.YMax)
     Zmid = 0.5*(self.BoundBox.ZMin+self.BoundBox.ZMax)
     LX   = (self.BoundBox.ZMin+self.BoundBox.XLength)
     LY   = (self.BoundBox.ZMin+self.BoundBox.YLength)
     LZ   = (self.BoundBox.ZMin+self.BoundBox.ZLength)
     PXMin = GEOUNED_Surface(('Plane',(FreeCAD.Vector(self.BoundBox.XMin,Ymid,Zmid),FreeCAD.Vector( 1,0,0),LY,LZ)),self.BoundBox)
     PXMax = GEOUNED_Surface(('Plane',(FreeCAD.Vector(self.BoundBox.XMax,Ymid,Zmid),FreeCAD.Vector(-1,0,0),LY,LZ)),self.BoundBox)
     PYMin = GEOUNED_Surface(('Plane',(FreeCAD.Vector(Xmid,self.BoundBox.YMin,Zmid),FreeCAD.Vector(0, 1,0),LZ,LX)),self.BoundBox)
     PYMax = GEOUNED_Surface(('Plane',(FreeCAD.Vector(Xmid,self.BoundBox.YMax,Zmid),FreeCAD.Vector(0,-1,0),LZ,LX)),self.BoundBox)
     PZMin = GEOUNED_Surface(('Plane',(FreeCAD.Vector(Xmid,Ymid,self.BoundBox.ZMin),FreeCAD.Vector(0,0, 1),LX,LY)),self.BoundBox)
     PZMax = GEOUNED_Surface(('Plane',(FreeCAD.Vector(Xmid,Ymid,self.BoundBox.ZMax),FreeCAD.Vector(0,0,-1),LX,LY)),self.BoundBox)

     return (PXMin,PXMax,PYMin,PYMax,PZMin,PZMax)

  def __removeExtraComp__(self,Obj,Box,mode='box'):
      reducedSol = []
      reducedDef = BoolSequence(operator='OR')
      if not Obj.Solids: return
      # Compare Solid BoundBox (here Box is BoundBox Object)
      if mode == 'box' :
        for i,sol in enumerate(Obj.Solids):
          if sol.BoundBox.isValid() :
            if Box.intersect(sol.BoundBox) :
              reducedSol.append(sol)
              reducedDef.append(Obj.Definition.elements[i])

      # Compare solid using distToshape (here Box is a Solid Cube object)       
      else :
        for i,sol in enumerate(Obj.Solids):
          dist = Box.distToShape(sol)[0]
          if dist == 0 :
            reducedSol.append(sol)
            reducedDef.append(Obj.Definition.elements[i])

      if len(reducedSol) < len(Obj.Solids):      
        Obj.updateSolids(reducedSol)
        Obj.setDefinition(reducedDef)
      return
            

  def __copyMeta__(self,m):
      solidsCopy = m.Solids[:]
      facesCopy = m.Faces[:]
      Meta = GEOUNED_Solid(m.__id__,solidsCopy) 
      Meta.setDefinition(m.Definition.copy())
      Meta.setFaces(facesCopy)
      if m.IsEnclosure :
          Meta.IsEnclosure = True
          Meta.EnclosureID = m.EnclosureID
      return Meta



    
