import os
import sys
import re
import FreeCAD
from GEOReverse.Modules.Parser import parser  as mp
from GEOReverse.Modules.remh import cell_card_string, remove_hash
from GEOReverse.Modules.Objects import *
import math
import numpy as np
from numpy import linalg as LA

class MCNPinput:
    def __init__(self,name):
        if not os.path.isfile(name):
            print("File %s does not exist" %name)
            sys.exit()
            return
        self.__inputcards__  = list(mp.get_cards(name))
        self.Transformations = self.__getTransList__()
        return

    def GetFilteredCells(self,Surfaces,config):
        levels,contLevels,Universes = self.GetLevelStructure()

        FilteredCells = {} 
        
        Ustart = config['Ustart'] 
        subUniverses = getSubUniverses(Ustart,Universes)
        subUniverses.add(Ustart)

        if config['levelMax'] =='all':
           levelMax = len(levels)
        else :   
           levelMax = config['levelMax'] + 1             

        levelUniverse = set() 
        for lev in range(0,levelMax) :             
           for U in levels[lev]:
              levelUniverse.add(U)
        subUniverses = subUniverses.intersection(levelUniverse)

        for U in list(Universes.keys()):
            if U not in subUniverses : del(Universes[U])
              
        for U in Universes.keys():
            FilteredCells[U] = selectCells(Universes[U],config)
            processSurfaces(FilteredCells[U],Surfaces)

        # change the surface name in surface dict
        newSurfaces = {}
        for k in Surfaces.keys():
           newkey = Surfaces[k].id
           newSurfaces[newkey] = Surfaces[k]

        
        for U,universe in FilteredCells.items() :
            substituteLikeCell(universe,newSurfaces)
            
            # set cell as CAD cell Object
            for cname,c in universe.items():
               #print(cname,c.geom.str)
               universe[cname] = CADCell(c)
           
        return levels,FilteredCells,newSurfaces     



    def GetLevelStructure(self):      
        containers = []
        Universe_dict = {}
        cellCards = {}

        
                
        for c in self.__inputcards__ :
          if c.ctype != mp.CID.cell : continue
          c.get_values()
          cstr = cell_card_string(''.join(c.lines))

          if cstr.TRCL : cstr.TRCL = TransformationMatrix(cstr.TRCL,self.Transformations)
          if cstr.TR   : cstr.TR   = TransformationMatrix(cstr.TR,self.Transformations)
          cellCards[c.name] = cstr
 
        setExplicitCellDefinition(cellCards)  
        for cname,c in cellCards.items():  
          if c.U is None: c.U = 0
          if c.U not in Universe_dict.keys(): Universe_dict[c.U] = {}
          Universe_dict[c.U].update({c.name:c})

          if c.FILL :
             containers.append(c)
             
        
        currentLevel = [0]
        nextLevel    = []
        contLevel    ={0:[(0,0)]}
        univLevel    = {0 : {0}}
        level = 0
        
        while True:
           level += 1 
           contLevel[level] = []
           univLevel[level] = set()
           for c in reversed(containers) :
              if c.U in currentLevel :
                 c.Level = level
                 nextLevel.append(c.FILL)
                 contLevel[level].append((c.U,c.name))
                 univLevel[level].add(c.FILL)
                 
                 containers.remove(c)

           if nextLevel == [] : break
           currentLevel = nextLevel
           nextLevel = []
 
        return univLevel,contLevel,Universe_dict   

    
    def GetCells(self,U=None,Fill=None):
        cell_cards = {}
        for c in self.__inputcards__ :
          if c.ctype != mp.CID.cell : continue
          c.get_values()
          U_cell = c.get_u()
          Fill_cell = c.get_f()
          if ( U is None and Fill is None) :
            cell_cards[c.name] = c
          elif U_cell == U and U is not None:
            cell_cards[c.name] = c  
          elif Fill_cell == Fill and Fill is not None:
            cell_cards[c.name] = c  

        return cell_cards  
        
    def GetSurfaces(self,scale=1.):
        surf_cards = {}
        number = 1
        for c in self.__inputcards__ :
          if c.ctype != mp.CID.surface : continue
          c.get_values()
          c.TR = TransformationMatrix(c.TR,self.Transformations)
          surf_cards[c.name] = (c.stype,c.scoefs,c.TR,number)
          number += 1

        # return surface as surface Objects type
        return Get_primitive_surfaces(surf_cards,scale)

    def __getTransList__(self):
        trl = {}
        for c in self.__inputcards__ :
          if c.ctype != mp.CID.data : continue
          c.get_values()
          if c.dtype == 'TRn' :
             trValues = []
             for v in c.values[1:]:
                trValues.append(v[0])
             trl[c.name]= getTransMatrix(trValues,c.unit)     
        return trl     


def getTransMatrix(trsf,unit='',scale=10.):
    
    if len(trsf) == 3 :
       trsfMat = FreeCAD.Matrix(1,0,0, trsf[0]*scale,
                                0,1,0, trsf[1]*scale,
                                0,0,1, trsf[2]*scale,
                                0,0,0, 1)
    else:
       if unit == '*':
          coeff = tuple(map(math.radians,trsf[3:12] )) 
          coeff = tuple(map(math.cos,coeff))
       else:
          coeff = trsf[3:12] 

       trsfMat = FreeCAD.Matrix(coeff[0],coeff[3] ,coeff[6], trsf[0]*scale,
                                coeff[1],coeff[4] ,coeff[7], trsf[1]*scale,
                                coeff[2],coeff[5] ,coeff[8], trsf[2]*scale,
                                0,0,0, 1)
    
def substituteLikeCell(universe,Surfaces):
    number = re.compile(r"\#?\s*\d+")
    newId = len(Surfaces)

    # create new cells objects for like cell cells    
    for c in universe.values(): 
       if c.likeCell : c.geom = universe[c.likeCell].geom.copy()
       if not c.TRCL : continue   # likebut cell should have TRCL card

    # transform change cell the parameters if needed
    for c in universe.values():
       if not c.TRCL : continue 
       cellSurf = c.geom.getSurfacesNumbers()
       surfDict = {}
       for surf in cellSurf:
           newId += 1
           surfDict[surf] = newId
           Surfaces[newId] = Surfaces[surf].copy()
           Surfaces[newId].id = newId
           Surfaces[newId].transform(c.TRCL)
       if c.FILL : c.TR = c.TRCL
       c.TRCL = None
           
       # rename the surface of the like cell
       pos = 0
       while True:
         m = number.search(c.geom.str,pos)
         if not m : break
         if "#" in m.group():
               pos = m.end()
               continue
         surf = int(m.group())
         pos = c.geom.replace(surf,surfDict[surf],pos)
         if pos < 0 : break            
       
def selectCells(cellList,config):
    selected = {}
    # options are 'all' material                     
    if config['mat'][0] == 'all' : 
       if config['cell'][0]=='all' :
          selected = cellList
       elif config['cell'][0]=='exclude' :
          for name,c in cellList.items(): 
             if name not in config['cell'][1]:
                selected[name] = c
       elif config['cell'][0]=='include' :
          for name,c in cellList.items(): 
             if name in config['cell'][1]:
                selected[name] = c

     # options are 'exclude' material
    elif config['mat'][0] == 'exclude':            
       if config['cell'][0]=='all' :
          for name,c in cellList.items():
             if c.FILL is None:
                 if c.MAT not in config['mat'][1]: 
                    selected[name] = c
             else:       
                 selected[name] = c # Fill cell are not tested against material number
       elif config['cell'][0]=='exclude' :
          for name,c in cellList.items():
             if c.FILL is None:
                 if c.MAT not in config['mat'][1]:
                    if name not in config['cell'][1]: 
                      selected[name] = c
             else:
                 if name not in config['cell'][1]:
                    selected[name] = c # Fill cell are not tested against material number
       elif config['cell'][0]=='include' :
          for name,c in cellList.items():
             if c.FILL is None:
                 if c.MAT not in config['mat'][1]:
                    if name in config['cell'][1]: 
                      selected[name] = c
             else:
                 if name in config['cell'][1]:
                    selected[name] = c # Fill cell are not tested against material number

     # options are 'include' material
    elif config['mat'][0] == 'include':            
       if config['cell'][0]=='all' :
          for name,c in cellList.items():
             if c.FILL is None:
                 if c.MAT in config['mat'][1]: 
                    selected[name] = c
             else:       
                 selected[name] = c # Fill cell are not tested against material number
       elif config['cell'][0]=='exclude' :
          for c in cellList:
             if c.FILL is None:
                 if c.MAT in config['mat'][1]:
                    if name not in config['cell'][1]: 
                      selected[name] = c
             else:
                 if name not in config['cell'][1]:
                    selected[name] = c # Fill cell are not tested against material number
       elif config['cell'][0]=='include' :
          for name,c in cellList.items():
             if c.FILL is None:
                 if c.MAT in config['mat'][1]:
                    if name in config['cell'][1]: 
                      selected[name] = c
             else:
                 if name in config['cell'][1]:
                    selected[name] = c # Fill cell are not tested against material number


    # remove complementary in cell of the universe
    for cname,c in selected.items() :
       c.geom = remove_hash(cellList,cname)

    return selected    
        
#Change implicit cell definition (like-but type or cell with #)
# to explicit cell definition (only surface number)
def setExplicitCellDefinition(universeCells):
    for cname,c in universeCells.items():
        if c.likeCell :
           lkc = universeCells[c.likeCell]
           c.geom = lkc.geom
           if not c.U     : c.U    = lkc.U
           if not c.FILL  : c.FILL = lkc.FILL
           if not c.MAT   : c.MAT  = lkc.MAT
           if not c.TR    : c.TR   = lkc.TR
           if not c.TRCL  : c.TRCL = lkc.TRCL
    return         


def processSurfaces(UCells,Surfaces):          
    number = re.compile(r"\#?\s*\d+")

    for cname,c in UCells.items() :
        c.geom.remove_comments(full=True) 
        pos = 0
        while True:
           m = number.search(c.geom.str,pos)
           if not m : break
           if "#" in m.group():
               pos = m.end()
               continue
           surf = int(m.group())
           if surf == 0 :
              print(c.name)
              print(m)
              print(c.geom.str)
           pos = c.geom.replace(surf,Surfaces[surf].id,pos)


def getTransMatrix(trsf,unit='',scale=10.):

    if len(trsf) == 3 :
       trsfMat = FreeCAD.Matrix(1,0,0, trsf[0]*scale,
                                0,1,0, trsf[1]*scale,
                                0,0,1, trsf[2]*scale,
                                0,0,0, 1)
    else:
       if unit == '*':
          coeff = tuple(map(math.radians,trsf[3:12] )) 
          coeff = tuple(map(math.cos,coeff))
       else:
          coeff = trsf[3:12] 

       trsfMat = FreeCAD.Matrix(coeff[0],coeff[3] ,coeff[6], trsf[0]*scale,
                                coeff[1],coeff[4] ,coeff[7], trsf[1]*scale,
                                coeff[2],coeff[5] ,coeff[8], trsf[2]*scale,
                                0,0,0, 1)
    return trsfMat   
       
def TransformationMatrix(TRSF,Transformations):

    if TRSF :
       if type(TRSF) is int :
          return Transformations[TRSF]

       else:
          return getTransMatrix(TRSF) 
    else:      
       return TRSF      



def getSubUniverses(Ustart,Universes):
    Uid = set() 
    for c in Universes[Ustart].values() :
        if c.FILL : Uid.add(c.FILL)

    AllU = Uid.copy()
    for U in Uid:
        AllU = AllU.union(getSubUniverses(U,Universes))
        
    return AllU    

# traduce mcnp surface definition for Solid_Cell class
#  planes:
#     Stype = 'plane'
#     params = [ax,ay,az,d]
#
#  spheres:
#     Stype = 'shpere'
#     params = [cx,cy,cz,R]
#
#  cylinders:
#     Stype = 'cylinder'
#     params = [[px,py,pz],[vx,vy,vz],R]
#
#  cones:
#     Stype = 'cone'
#     params = [[px,py,pz],[vx,vy,vz],t,sht]
#
#  torus:
#     Stype = 'torus'
#     params = [[px,py,pz],[vx,vy,vz],ra,r]

# Return a diccionary with the corresponding surface Object
def Get_primitive_surfaces(mcnp_surfaces,scale=10.) :

      X_vec = FreeCAD.Vector(1.0,0.0,0.0)
      Y_vec = FreeCAD.Vector(0.0,1.0,0.0)
      Z_vec = FreeCAD.Vector(0.0,0.0,1.0)
      negX_vec = -X_vec
      negY_vec = -Y_vec
      negZ_vec = -Z_vec
      origin= FreeCAD.Vector(0.0,0.0,0.0)

      surfaces = {}
      for Sid in mcnp_surfaces.keys() :      
         MCNPtype   = mcnp_surfaces[Sid][0].upper()
         MCNPparams = mcnp_surfaces[Sid][1]
         trsf       = mcnp_surfaces[Sid][2]   
         number     = mcnp_surfaces[Sid][3]
         
         params = []
         Stype = None
         if MCNPtype in ('P','PX','PY','PZ'):
            Stype = 'plane'
            if MCNPtype == 'P':
                if len(MCNPparams) == 4 :
                   normal = FreeCAD.Vector(MCNPparams[0:3])
                   params = (normal, MCNPparams[3]*scale)
                else:
                   coeffs = pointsToCoeffs(MCNPparams[0:9])
                   normal = FreeCAD.Vector(coeffs[0:3])
                   point  = coeffs[3]/normal.Length
                   normal.normalize()
                   params = (normal, point*scale)
            elif MCNPtype == 'PX':
                params = (X_vec, MCNPparams[0]*scale )
            elif MCNPtype == 'PY':
                params = ( Y_vec, MCNPparams[0]*scale )
            elif MCNPtype == 'PZ':
                params = ( Z_vec, MCNPparams[0]*scale )
                
         elif MCNPtype in ['S','SX','SY','SZ','SO'] :
            Stype = 'sphere'
            if MCNPtype == 'S':
                params = (FreeCAD.Vector(MCNPparams[0:3])*scale, MCNPparams[3]*scale)
            elif MCNPtype == 'SX':
                params = ( FreeCAD.Vector(MCNPparams[0]*scale, 0.0, 0.0), MCNPparams[1]*scale )
            elif MCNPtype == 'SY':
                params = ( FreeCAD.Vector(0.0, MCNPparams[0]*scale, 0.0), MCNPparams[1]*scale )
            elif MCNPtype == 'SZ':
                params = ( FreeCAD.Vector(0.0, 0.0, MCNPparams[0]*scale), MCNPparams[1]*scale )
            elif MCNPtype == 'SO':
                params = ( origin, MCNPparams[0]*scale )

         elif MCNPtype in ['CX','C/X','CY','C/Y','CZ','C/Z'] :
            if MCNPtype in ['CX','CY','CZ']:
               R = MCNPparams[0]
            else:
               R = MCNPparams[2]
               x1 = MCNPparams[0]
               x2 = MCNPparams[1]
               
            Stype = 'cylinder'
        
            if MCNPtype == 'CX':
                v = X_vec
                p = origin
            elif MCNPtype == 'CY':
                v = Y_vec
                p = origin
            elif MCNPtype == 'CZ':
                v = Z_vec
                p = origin
            elif MCNPtype ==  'C/X':
                v = X_vec
                p = FreeCAD.Vector(0., x1, x2)
            elif MCNPtype == 'C/Y':
                v = Y_vec
                p = FreeCAD.Vector(x1, 0., x2)
            elif MCNPtype == 'C/Z':
                v = Z_vec
                p = FreeCAD.Vector(x1, x2, 0.)

            if scale != 1.0 :
              p = p.multiply(scale)
              R *= scale

            params = ( p,v,R )  

         elif MCNPtype in ['KX','KY','KZ'] :
            Stype = 'cone' 
            x1 = MCNPparams[0]
            t2 = MCNPparams[1]
            t  = math.sqrt(t2)
            dblsht = True

            if len(MCNPparams) == 3:
              dblsht = False
              sht = MCNPparams[2]  
              if sht == 0. :
                dblsht = True

            if MCNPtype == 'KX':
                p  = FreeCAD.Vector(x1,0.0,0.0)
                v  = X_vec
                if not dblsht :
                  if sht < 0 :
                    v = negX_vec
            elif MCNPtype == 'KY':
                p  = FreeCAD.Vector(0.0,x1,0.0)
                v  = Y_vec
                if not dblsht :
                  if sht < 0 :
                     v = negY_vec
            elif MCNPtype == 'KZ':
                p  = FreeCAD.Vector(0.0,0.0,x1)
                v  = Z_vec
                if not dblsht :
                  if sht < 0 :
                     v = negZ_vec             

            p = p.multiply(scale)
            params = (p,v,t,dblsht)

                         
         elif MCNPtype in ['K/X','K/Y','K/Z'] :
            Stype = 'cone' 
            x1 = MCNPparams[0]
            x2 = MCNPparams[1]
            x3 = MCNPparams[2]
            p  = FreeCAD.Vector(x1,x2,x3)
            t2 = MCNPparams[3]
            t  = math.sqrt(t2)
            dblsht = True

            if len(MCNPparams) == 5:
              dblsht = False
              sht = MCNPparams[4]  
              if sht == 0. :
                dblsht = True

            if MCNPtype == 'K/X':
              v  = X_vec
              if not dblsht :
                 if sht < 0 :
                    v = negX_vec
            elif MCNPtype == 'K/Y':
                v  = Y_vec
                if not dblsht :
                   if sht < 0 :
                      v = negY_vec
            elif MCNPtype == 'K/Z':
                v  = Z_vec
                if not dblsht :
                   if sht < 0 :
                      v = negZ_vec             

            p = p.multiply(scale)
            params = (p,v,t,dblsht)

         elif MCNPtype in ['TX','TY','TZ'] :
            Stype = 'torus'
            x1 = MCNPparams[0]
            x2 = MCNPparams[1]
            x3 = MCNPparams[2]
            p  = FreeCAD.Vector(x1,x2,x3)
            Ra = MCNPparams[3]
            r1 = MCNPparams[4]
            r2 = MCNPparams[5]
            if (r1 != r2 ) :
               print('ellipsoid torus not implemented : {} {}'.format(r1,r2) )
            R = (r1+r2)*0.5         

            if   MCNPtype == 'TX':
              v  = X_vec
            elif MCNPtype == 'TY':
              v  = Y_vec
            elif MCNPtype == 'TZ':
              v  = Z_vec
                
            if scale != 1.0 :
              Ra *= scale
              R  *= scale
              p = p.multiply(scale)

            params = ( p,v,Ra,R)  

         elif MCNPtype == 'GQ' :
            Qparams = tuple(MCNPparams[0:10])   
            Stype,quadric = gq2cyl(Qparams)

            if Stype == 'cylinder' :
                p = FreeCAD.Vector(quadric[0:3])
                v = FreeCAD.Vector(quadric[3:6])
                R = quadric[6]
                if scale != 1.0 :
                  R *= scale
                  p = p.multiply(scale)

                params = ( p,v,R )
                
            elif Stype == 'cone' :
                p = FreeCAD.Vector(quadric[0:3])
                v = FreeCAD.Vector(quadric[3:6])
                t = quadric[6]
                dblsht = quadric[7]
                if scale != 1.0 :
                    p = p.multiply(scale)
                 
                params = ( p,v,t,dblsht )

            else:
                print ( Stype )
                params = None
#                get_quadric_surface(params)

         elif MCNPtype == 'X' :
             if len(MCNPparams) == 2 :
                Stype = 'plane'
                params = (X_vec, MCNPparams[0]*scale )
             elif len(MCNPparams) == 4 :
                 if ( abs(MCNPparams[1] - MCNPparams[3]) ) > 1.e-12 :
                    Stype = 'cone'
                    dblsht = False
                    t = (MCNPparams[3] - MCNPparams[1])/(MCNPparams[2] - MCNPparams[0])
                    x =  MCNPparams[0] - MCNPparams[1]/t 
                    if (MCNPparams[0]-x)*(MCNPparams[2]-x) > 0 :
                        p = FreeCAD.Vector(x,0.,0.)
                        if (MCNPparams[0]-x) > 0 :
                           v = X_vec
                        else:
                           v = negX_vec 
                        if scale != 1.0 : p *= scale
                        params = (p,v,abs(t),dblsht)
                 elif  abs(MCNPparams[1]) < 1.e-12 :
                    Stype = 'plane' 
                    if ( abs(MCNPparams[0] - MCNPparams[2]) ) < 1.e-12 :
                        params = (X_vec, MCNPparams[0]*scale )
                 else:    
                    Stype = 'cylinder'
                    if scale != 1.0 :
                       p = p.multiply(scale)
                       R *= scale
                    params = ( origin,X_vec,MCNPparams[1] )
             else:
                 print('not implemented surfaces defined by point with more than 2couples of value')  

         elif MCNPtype == 'Y' :
             if len(MCNPparams) == 2 :
                Stype = 'plane'
                params = (Y_vec, MCNPparams[0]*scale )
             elif len(MCNPparams) == 4 :
                 if ( abs(MCNPparams[1] - MCNPparams[3]) ) > 1.e-12 :
                    Stype = 'cone'
                    dblsht = False
                    t = (MCNPparams[3] - MCNPparams[1])/(MCNPparams[2] - MCNPparams[0])
                    y = MCNPparams[0] - MCNPparams[1]/t
                    if (MCNPparams[0]-y)*(MCNPparams[2]-y) > 0 :
                        p = FreeCAD.Vector(0.,y,0.)
                        if (MCNPparams[0]-y) > 0 :
                           v = Y_vec
                        else:
                           v = negY_vec 
                        if scale != 1.0 : p = p.multiply(scale)
                        params = (p,v,abs(t),dblsht)
                 elif  abs(MCNPparams[1]) < 1.e-12 :
                    Stype = 'plane' 
                    if ( abs(MCNPparams[0] - MCNPparams[2]) ) < 1.e-12 :
                        params = (Y_vec, MCNPparams[0]*scale )
                 else:    
                    Stype = 'cylinder'
                    if scale != 1.0 :
                       p = p.multiply(scale)
                       R *= scale
                    params = ( origin,Y_vec,MCNPparams[1] )
             else:
                 print('not implemented surfaces defined by point with more than 2couples of value')
                 
         elif MCNPtype == 'Z' :
             if len(MCNPparams) == 2 :
                Stype = 'plane'
                params = (Z_vec, MCNPparams[0]*scale )
             elif len(MCNPparams) == 4 :
                 if ( abs(MCNPparams[1] - MCNPparams[3]) ) > 1.e-12 :
                    Stype = 'cone'
                    dblsht = False
                    t = (MCNPparams[3] - MCNPparams[1])/(MCNPparams[2] - MCNPparams[0])
                    z = MCNPparams[0] - MCNPparams[1]/t
                    if (MCNPparams[0]-z)*(MCNPparams[2]-z) > 0 :
                        p = FreeCAD.Vector(0.,0.,z)
                        if (MCNPparams[0]-z) > 0 :
                           v = Z_vec
                        else:
                           v = negZ_vec 
                        if scale != 1.0 : p = p.multiply(scale)
                        params = (p,v,abs(t),dblsht)
                 elif  abs(MCNPparams[1]) < 1.e-12 :
                    Stype = 'plane' 
                    if ( abs(MCNPparams[0] - MCNPparams[2]) ) < 1.e-12 :
                        params = (Z_vec, MCNPparams[0]*scale )
                 else:    
                    Stype = 'cylinder'
                    if scale != 1.0 :
                       p = p.multiply(scale)
                       R *= scale
                    params = ( origin,Z_vec,MCNPparams[1] )
             else:
                 print('not implemented surfaces defined by point with more than 2couples of value')                 
     

         if Stype == 'plane':
             surfaces[Sid] = Plane(number, params, trsf)
         elif Stype == 'sphere':
             surfaces[Sid] = Sphere(number, params, trsf)    
         elif Stype == 'cylinder':
             surfaces[Sid] = Cylinder(number, params,trsf)
         elif Stype == 'cone':
             surfaces[Sid] = Cone(number, params,trsf)
         elif Stype == 'torus':
             surfaces[Sid] = Torus(number, params,trsf)
         else :
             print('Undefined',Sid)
             print( MCNPtype ,number,MCNPparams) 

      return surfaces 

def gq2cyl(x):
# Conversion de GQ a Cyl
# Ax2+By2+Cz2+Dxy+Eyz+Fxz+Gx+Hy+Jz+K=0
# x.T*M*x + b.T*x + K = 0
  minWTol=5.e-2
  minRTol=1.e-3
  #minRTol=3.e-1
  # lx = np.array(x)
  tp = ''
  M = np.array( [[x[0],x[3]/2,x[5]/2], \
                 [x[3]/2,x[1],x[4]/2], \
                 [x[5]/2,x[4]/2,x[2]]] )
  w,P = LA.eigh(M)
  sw = np.sort(w)
  aw = np.abs(w)
  asw= np.sort(aw)
  # Test for cylinder (least abs value is much less than others)
  if asw[0]<minWTol*asw[1]:
    tp = 'cylinder'
    rv = [0.]*7   # X0,Y0,Z0, VX, VY, VZ, R
    iaxis = np.where(aw==asw[0])[0][0]
    otherAxes = ( (iaxis+1)%3, (iaxis+2)%3 )
    if abs(w[otherAxes[0]]-w[otherAxes[1]])>minRTol*asw[2]:
       tp = 'not found - ellipsoid cylinder'
       rv = [0]
       return tp,rv
    # Vector de desplazamiento
    # x0 = -0.5*Pt*D-1*P*b pero ojo que un lambda es cero
    # P es la matriz de valores propios
    b = np.array(x[6:9])
    Pb= np.matmul(P,b)
    for i in otherAxes: Pb[i] /= w[i]
    x0 = -0.5*np.matmul(P.T,Pb)
    k  = -0.5*np.matmul(x0,b) - x[9]
    # Resultados finales

    rv[0:3] = x0                     # Punto del eje
    rv[3:6] = P[:,iaxis]             # Vector director
    rv[6]   = np.sqrt(k/sw[1])       # Radio
  # Test for cone (incomplete, returns empty data list)
  elif np.sign(sw[0])!=np.sign(sw[2]):   # maybe cone
    tp = 'cone'
    rv = [0.]*8   #  X0, Y0, Z0, VX, VY, VZ, tgAlpha, double sheet
    if np.sign(sw[0])==np.sign(sw[1]):
      iaxis = np.where(w==sw[2])[0][0]
    else:
      iaxis = np.where(w==sw[0])[0][0]
    otherAxes = ( (iaxis+1)%3, (iaxis+2)%3 )
    if abs(w[otherAxes[0]]-w[otherAxes[1]])>minRTol*asw[2]:
       tp = 'not found - ellipsoid cone/hyperboloid'
       rv = [0]
       return tp,rv
    # Displacement vector ( x0 = -0.5*M^-1*b = -0.5*P.T*D^-1*P*b
    b = np.array(x[6:9])
    x0= -0.5*np.matmul(P,np.matmul(P.T,b)/w)
    k = x0.T @ M @ x0 - x[9]  
    if np.abs(k*w[iaxis])>minRTol*minRTol*asw[2]:
       
      # tp = 'not found - hyperboloid'
      # rv = [0]
        # force cone surface
        print('Force cone surface')
        tp = 'cone'
        rv[0:3] = x0                                   # vertex point
        rv[3:6] = P[:,iaxis]                           # axis direction
        rv[6]   = np.sqrt(-w[otherAxes[0]]/w[iaxis])   # semiangle tangent
        rv[7]   = True                                 # here always double sheet cones
        return tp,rv   
    # return value
    rv[0:3] = x0                                   # vertex point
    rv[3:6] = P[:,iaxis]                           # axis direction
    rv[6]   = np.sqrt(-w[otherAxes[0]]/w[iaxis])   # semiangle tangent
    rv[7]   = True                                 # here always double sheet cones
  else:
       tp = 'not found - unknown'
       rv = [0]
  return tp,rv

def pointsToCoeffs(scf):
    # mcnp implementation to convert 3 point plane to
    # plane parameters

    tpp = [0]*4
    for i in range(1,4) :
       j = i%3 + 1
       k = 6 -i -j
       k -= 1
       j -= 1
       tpp[i-1] =  scf[j  ]*(scf[k+3]-scf[k+6]) \
                 +scf[j+3]*(scf[k+6]-scf[k  ]) \
                 +scf[j+6]*(scf[k  ]-scf[k+3])
       tpp[3]  += scf[i-1]*(scf[j+3]*scf[k+6]-scf[j+6]*scf[k+3])

    xm = 0
    coeff = [0]*4
    for i in range(1,5):
       if  xm == 0 and tpp[4-i] != 0  : xm = 1/tpp[4-i]
       coeff[4-i] = tpp[4-i]*xm

    # coeff [0:3] a,b,c plane parameters
    # coeff [3]   d plane parameter
    # normalization is d set to one if origin is not in the plane
    return coeff

