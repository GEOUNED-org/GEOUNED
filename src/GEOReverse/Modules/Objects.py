from GEOReverse.Modules.buildSolidCell import BuildSolid
import FreeCAD, Part
import numpy as np
from GEOReverse.Modules.remh import cline
from GEOReverse.Modules.Utils.booleanFunction import outterTerms,BoolSequence
from GEOReverse.Modules.Utils.BooleanSolids import buildCTableFromSolids

class CADCell:
    def __init__(self,stringCell=None):

        if not stringCell :
            self.surfaces = {}
            self.surfaceList = []
            self.shape    =  None
            #self.likeCell = None
            self.definition = None
            self.name     = 0
            #self.TRCL     = None  # cell transformacion "like-but" cells
            self.TRFL     = None    # Universe transformation in fill Universe
            self.U        = -1     # Cell Universe number
            self.FILL     = 0  # Fill Universe number
            self.MAT      = 0   # material number
            self.CurrentTR = None
            self.level   = None
            self.__defTerms__ = None
            self.__operator__ = None
        else:
            self.surfaces = None
            self.shape    = None
            self.name     = stringCell.name
            self.TRFL     = stringCell.TR    # Universe transformation in fill Universe
            self.U        = stringCell.U     # Cell Universe number
            self.FILL     = stringCell.FILL  # Fill Universe number
            self.MAT      = stringCell.MAT   # material number
            self.CurrentTR = self.TRFL
            self.level  = None
            
            self.__defTerms__ = None
            self.__operator__ = None        
            self.__setDefinition__(stringCell)
               
    def copy(self):
        cpCell = CADCell() 
        cpCell.surfaceList = self.surfaceList[:]
        cpCell.surfaces = {}
        for name,s in self.surfaces.items():
            cpCell.surfaces[name] = s.copy()
            
        if type(self.definition) is cline :  
           cpCell.definition = cline(self.definition.str)

        elif type(self.definition) is BoolSequence:
           cpCell.definition = self.definition.copy()
         
        cpCell.name       = self.name
        cpCell.TRFL       = self.TRFL
        cpCell.U          = self.U
        cpCell.FILL       = self.FILL
        cpCell.MAT        = self.MAT
        cpCell.level      = self.level

        if self.CurrentTR is not None :
            cpCell.CurrentTR  = self.CurrentTR.submatrix(4)

        if self.shape is not None :
          cpCell.shape      = self.shape.copy()

        return cpCell


    def getSubCell(self,seq):

        subCell = self.copy()
        subCell.definition = seq.copy()
        subCell.shape = None
        subCell.surfaceList  = subCell.definition.getSurfacesNumbers()
        for s in tuple(subCell.surfaces.keys()) :
            if s not in subCell.surfaceList: del(subCell.surfaces[s])
                  
        return subCell    

    def split(self,nparts=2):

        if nparts == 1:
            return (self,None)
        terms,operador = self.getOuterTerms()
        nelemts = int(len(terms)/nparts)
        subDefList = []

        if operador == 'AND':
          for i in range(parts-1):
             newdef = ') ('.join(terms[i*nelemts:(i+1)*nelemts])
             newdef = '({})'.format(newdef)
             subDefList.append(newdef) 
          newdef = ') ('.join(terms[(nparts-1)*nelemts:])
          newdef = '({})'.format(newdef)
          subDefList.append(newdef)

        else:
          for i in range(nparts-1):
             newdef = '):('.join(terms[i*nelemts:(i+1)*nelemts])
             newdef = '({})'.format(newdef)
             subDefList.append(newdef)
          newdef = '):('.join(terms[(nparts-1)*nelemts:])
          newdef = '({})'.format(newdef)
          subDefList.append(newdef)

   
        subCellList=[] 
        for df in subDefList:
           subCell = self.copy()
           subCell.definition= cline(df)
           subCell.shape = None
           subCell.surfaceList  = subCell.definition.getSurfacesNumbers()
           for s in tuple(subCell.surfaces.keys()) :
               if s not in subCell.surfaceList: del(subCell.surfaces[s])
                  
           subCellList.append(subCell)    

        return subCellList,operador

    def getOuterTerms(self):
        if not self.__defTerms__ :
           self.__defTerms__,self.__operator__ = outterTerms(self.definition.str)
        return self.__defTerms__, self.__operator__
        
    def makeBox(self,boundBox):
        box_origin = FreeCAD.Vector(boundBox.XMin,boundBox.YMin,boundBox.ZMin)
        return Part.makeBox(boundBox.XLength,boundBox.YLength,boundBox.ZLength,box_origin )

    def buildShape(self,boundBox,force=False,surfTR=None,simplify=False,fuse=False):

        if self.shape is not None and not force: return
        if surfTR                  : self.transformSurfaces(surfTR)
        
        cutShape = BuildSolid(self,boundBox,simplify=simplify)
        
        if fuse or True:
          self.shape = FuseSolid(cutShape)
        else:
          self.shape = Part.makeCompound(cutShape)

    def buildSurfaceShape(self,boundBox):
        for s in self.surfaces.values():
           s.buildShape(boundBox)
        
    def transformSolid(self,matrix,reverse=False):
       if not self.shape : return 
       if reverse :
          self.shape = self.shape.transformGeometry(matrix.inverse())
       else:
          self.shape = self.shape.transformGeometry(matrix)

    def transformSurfaces(self,matrix):
        for s in self.surfaces.values():
            s.transform(matrix)
        
    def setSurfaces(self,Surfaces):
        if self.surfaces is not None : return
        self.surfaces = {}
        for s in self.surfaceList :
            self.surfaces[s] = Surfaces[s]

    def cleanUndefined(self):
       undefined = []
       for s in self.definition.getSurfacesNumbers() :
          if self.surfaces[s].params is None:
              undefined.append(s)
       if undefined:
           print('undefined surfaces',undefined)
           self.definition.removeSurface(undefined)

       for s in undefined:
           del (self.surfaces[s])
      
    def __setDefinition__(self,stringCell) :

       self.definition = stringCell.geom
       self.definition.remove_comments(full=True)
       self.definition.remove_cr()
       self.definition.remove_multispace()
       self.definition.remove_redundant()
       self.surfaceList  = self.definition.getSurfacesNumbers()
     

class Plane:
    def __init__(self,Id,params,tr=None):
        self.type = 'plane'
        self.id = Id
        self.shape = None
        self.params = params
        if tr :
           self.transform(tr)

    def __str__(self):
        return 'plane : {}\nParameters : {}'.format(self.id,self.params)

    def copy(self):
        return Plane(self.id,self.params)
    
    def transform(self,matrix):
        v,d = self.params
        p = d*v  # vector p is d*plane normal
        v = matrix.submatrix(3).multVec(v)
        v.normalize()
        d = matrix.multVec(p)*v
        self.params = (v,d)

    def buildShape(self,boundBox):
        normal, p0 = self.params
        Box = FreeCAD.BoundBox(boundBox)
        Box.enlarge(10)

        pointEdge=[]
        for i in range(12):
           edge = Box.getEdge(i)
           p1 = normal.dot(edge[0])
           p2 = normal.dot(edge[1])
           d0 = p0 - p1
           d1 = p2 - p1
           if d1 != 0 :
             a = d0/d1
             if (a >= 0 and a <= 1):
                pointEdge.append(edge[0]+a*(edge[1]-edge[0]))
  
        if len(pointEdge) == 0 : return
        s = FreeCAD.Vector((0,0,0))
        for v in pointEdge :
           s = s + v
        s = s /  len(pointEdge) 
   
        vtxvec = []
        for v in pointEdge :
           vtxvec.append(v-s)

        X0 = vtxvec[0]
        Y0 = normal.cross(X0)

        orden=[]
        for i,v in enumerate(vtxvec):
          phi = np.arctan2(v.dot(Y0),v.dot(X0))
          orden.append((phi,i))
        orden.sort()
    
        self.shape = Part.Face( Part.makePolygon([ pointEdge[p[1]] for p in orden ], True))
        

class Sphere:
    def __init__(self,Id,params,tr=None):
        self.type = 'sphere'
        self.id = Id
        self.shape = None
        self.params = params
        if tr :
           self.transform(tr)

    def copy(self):
        return Sphere(self.id,self.params)
    
    def transform(self,matrix):
        p,R = self.params
        p = matrix.multVec(p)
        self.params = (p,R)    

    def buildShape(self,boundBox):
        origin,R = self.params
        self.shape = Part.makeSphere(R,origin)
     

class Cylinder:
    def __init__(self,Id,params,tr=None):
        self.type = 'cylinder'
        self.id = Id
        self.shape = None
        self.params = params
        if tr :
           self.transform(tr)

    def copy(self):
        return Cylinder(self.id,self.params)

    def transform(self,matrix):
        p,v,R = self.params
        v = matrix.submatrix(3).multVec(v)
        p = matrix.multVec(p)
        self.params = (p,v,R)
    
    def buildShape(self,boundBox):

        p, vec, r = self.params

        dmin = vec.dot(boundBox.getPoint(0)-p)
        dmax = dmin
        for i in range(1,8) :
           d = vec.dot(boundBox.getPoint(i)-p)
           dmin = min(d,dmin)
           dmax = max(d,dmax)

        height = dmax-dmin
        dmin -= 0.1*height
        dmax += 0.1*height
        height = dmax-dmin
   
        point = p + dmin * vec
        self.shape = Part.makeCylinder( r,height,point,vec,360)
        return

class Cone:
    def __init__(self,Id,params,tr=None):
        self.type = 'cone'
        self.id = Id
        self.shape = None
        self.params = params
        if tr :
           self.transform(tr) 

    def copy(self):
        return Cone(self.id,self.params)

    def transform(self,matrix):
        p,v,t,dbl = self.params
        v = matrix.submatrix(3).multVec(v)
        p = matrix.multVec(p)
        self.params = (p,v,t,dbl)
    
    def buildShape(self,boundBox):
        apex, axis, t, dblsht = self.params
        
        dmin = axis.dot(boundBox.getPoint(0)-apex)
        dmax = dmin
        for i in range(1,8) :
           d = axis.dot(boundBox.getPoint(i)-apex)
           dmin = min(d,dmin)
           dmax = max(d,dmax)

        length = max(abs(dmin),abs(dmax))
        R = length * t
        OneSheetCone = Part.makeCone( 0,R,length,apex,axis,360)
        if not dblsht :
           self.shape = OneSheetCone
        else:
           OtherSheet = Part.makeCone( 0,R,length,apex,-axis,360)
           DoubleSheetCone = OneSheetCone.fuse([OtherSheet])
           DoubleSheetCone.removeSplitter()
           self.shape = DoubleSheetCone  

class Torus:
    def __init__(self,Id,params,tr=None):
        self.type = 'torus'
        self.id = Id
        self.shape = None
        self.params = params
        if tr :
           self.transform(tr)
           
    def copy(self):
        return Torus(self.id,self.params)

    def transform(self,matrix):
        p,v,Ra,R = self.params
        v = matrix.submatrix(3).multVec(v)
        p = matrix.multVec(p)
        self.params = (p,v,Ra,R)
    
    def buildShape(self,boundBox):
        center, axis, majorR, minorR  = self.params  #Ra distance from torus axis; R radius of toroidal-cylinder
        self.shape = Part.makeTorus(majorR,minorR,center,axis)
      

class Undefined:
    def __init__(self,Id):
        self.type = 'Undefined'
        self.id = Id
        self.shape = None
        self.params = None

    def copy(self):
        return Undefined(self.id)

    def buildShape(self,boundBox):
        return
    
    def transform(self,matrix):
        return

def FuseSolid(parts):
    if (len(parts)) <= 1:
       if parts : 
          solid = parts[0]
       else:
          return None 
    else:
       try:
          fused = parts[0].fuse(parts[1:])
       except:
          fused = None

       if fused is not None:
         try : 
             refinedfused = fused.removeSplitter()
         except :
             refinedfused = fused
           
         if refinedfused.isValid() :
             solid = refinedfused
         else :
             if fused.isValid():
                solid = fused
             else:
                solid = Part.makeCompound(parts)
       else:
         print('parts',parts)
         solid = Part.makeCompound(parts)
           
    if solid.Volume < 0 : solid.reverse()
    return solid
