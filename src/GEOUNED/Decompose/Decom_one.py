#   Conversion to MCNP v0.0 
#   Only one solid and planar surfaces
#

import FreeCAD, Part
import BOPTools.SplitAPI
import GEOUNED.Utils.Functions as UF
import GEOUNED.Utils.Geometry_GU as GU
import GEOUNED.Conversion.CellDefinition as CD
from collections import OrderedDict
from GEOUNED.Utils.BasicFunctions_part1 import isParallel, isSameValue,\
     isInLine, isInPlane
from GEOUNED.Utils.BasicFunctions_part2 import isDuplicateInList
from GEOUNED.Utils.Options.Classes import Tolerances as tol
from GEOUNED.Utils.Options.Classes import Options as opt
import math
import copy

twoPi = math.pi*2
    
def splitFullCylinder(solid):
    explode = []
    Bases = [solid]
    while True:
        newBases = []
        for base in Bases:
            cutSolids = cutFullCylinder(base)
            if len(cutSolids) == 1:
                explode.extend(cutSolids)
            else:
                newBases.extend(cutSolids)
        if len(newBases) == 0 :
            break
        else:
            Bases = newBases
         
    return Part.makeCompound(explode)

def cutFullCylinder(solid):
    solid_GU = GU.solid_GU(solid)
    Surfaces=UF.Surfaces_dict()    
    flag_inv = CD.isInverted(solid_GU.solid)
    UniverseBox = solid.BoundBox

    for face in solid_GU.Faces:
       surf = str(face.Surface)
       if surf == '<Cylinder object>' :
          if (flag_inv):
              orient = 'Reversed' if face.Orientation =='Forward' else 'Forward'
          else:
              orient=face.Orientation
    
          u1,u2,v1,v2 = face.ParameterRange
          angle = abs(u2-u1)

          # closed convex cylinder
          if abs(angle%twoPi) < 1e-2 and orient=='Forward':
             dir  = face.Surface.Axis
             orig = face.Surface.Center
             rad  = face.Surface.Radius
             dimL = face.ParameterRange[3]-face.ParameterRange[2]
             cylinder = UF.GEOUNED_Surface(('Cylinder',(orig,dir,rad,dimL)),UniverseBox)
             cylinder.buildSurface()
             Surfaces.addCylinder(cylinder,False)

             #add planes if cylinder axis is cut by a plane (plane quasi perpedicular to axis)
             for p in CylBoundPlanes(face,UniverseBox) :
                p.buildSurface()
                Surfaces.addPlane(p,False)
             break
                
    Planes = [] 
    for P in ('PX','PY','PZ','P') : Planes.extend(Surfaces[P])
    Planes = sortPlanes(Planes,True)

    
    if len(Planes) == 0 : return [solid]
    if len(Planes[-1]) < opt.nPlaneReverse : Planes.reverse()

    cut = False
    for pp in Planes :
       # cut with more external parallel planes
       pp[0].buildSurface()
       pp[-1].buildSurface()
       tools = (pp[0].shape, pp[-1].shape) 
       comsolid=BOPTools.SplitAPI.slice(solid,tools,'Split', tolerance = opt.splitTolerance)
       if len(comsolid.Solids) > 1 :
         outSolid = comsolid.Solids
         cut = True
         break

    if cut : return outSolid

    tool = (Surfaces['Cyl'][0].shape,)
    comsolid=BOPTools.SplitAPI.slice(solid,tool,'Split', tolerance = opt.splitTolerance)

    if len(comsolid.Solids) > 1 :
       outSolid = comsolid.Solids
    else:
       outSolid = [solid]

    return outSolid   


def GenPlane(pos,normal,diag):
    plane = Part.makePlane(diag,diag,pos,normal)
    vec_on_plane= plane.Vertexes[3].Point.sub(plane.Vertexes[0].Point)
    new_pos = plane.Vertexes[0].Point.sub(vec_on_plane)
    plane_center = Part.makePlane(2.0*diag,2.0*diag,new_pos,normal)
    return plane_center

def CylBoundPlanes(face,boundBox):
    Edges=face.OuterWire.Edges
    planes = []
    for e in Edges:
      try:
         curve = str(e.Curve)
      except:
         curve = 'none'
         
      if curve[0:6] == 'Circle' :
          dir = e.Curve.Axis
          center = e.Curve.Center
          dim1 = e.Curve.Radius
          dim2 = e.Curve.Radius
          plane = UF.GEOUNED_Surface(('Plane',(center,dir,dim1,dim2)),boundBox)
          planes.append(plane)
        
      elif curve == '<Ellipse object>' :    
          dir = e.Curve.Axis
          center = e.Curve.Center
          dim1 = e.Curve.MinorRadius
          dim2 = e.Curve.MajorRadius
          plane = UF.GEOUNED_Surface(('Plane',(center,dir,dim1,dim2)),boundBox)
          planes.append(plane)   

    return planes


def TorusBoundPlanes(face,boundBox):
    params = face.ParameterRange
    planes = []
    if isSameValue(params[1]-params[0],twoPi,tol.value) : return planes
    
    Edges=face.OuterWire.Edges
    
    for e in Edges:
      try:
         curve = str(e.Curve)
      except:
         curve = 'none'
   
      if curve[0:6] == 'Circle' :
          dir = e.Curve.Axis
          if not isParallel (dir,face.Surface.Axis,tol.angle) :
            center = e.Curve.Center
            dim1 = e.Curve.Radius
            dim2 = e.Curve.Radius
            plane = UF.GEOUNED_Surface(('Plane',(center,dir,dim1,dim2)),boundBox)
            planes.append(plane)
        
      elif curve == '<Ellipse object>' :    
          dir = e.Curve.Axis
          center = e.Curve.Center
          dim1 = e.Curve.MinorRadius
          dim2 = e.Curve.MajorRadius
          plane = UF.GEOUNED_Surface(('Plane',(center,dir,dim1,dim2)),boundBox)
          planes.append(plane)

      elif curve == '<BSplineCurve object>' :
          planeParams = PlaneSplineCurve(e)
          if planeParams is not None:
             plane = UF.GEOUNED_Surface(('Plane',planeParams),boundBox)
             planes.append(plane) 

    return planes

def PlaneSplineCurve(edge):

    normal = edge.derivative1At(0).cross(edge.derivative1At(0.5))
    normal.normalize()
    curve2D = True
    for p in (0.25,0.75,1) :
       # check if derivative orthogonal to curve normal vector 
       if abs(normal.dot(edge.derivative1At(p))) > tol.value :
           curve2D = False
           break

    r = edge.valueAt(0.25)-edge.valueAt(0.75)
    if curve2D:    
       return (edge.valueAt(0),normal,r.Length,r.Length)
    else:
       return None 


def ExtractSurfaces(solid,kind,UniverseBox,MakeObj=False):
    if kind == 'All':
       fuzzy = True 
       solidParts = []
       for s in solid.Solids:
         solidParts.append(GU.solid_GU(s))
    else:
        fuzzy = False
        if kind == 'Plane3Pts' :
           P3P = True
        else:
           P3P = False
        solidParts = [GU.solid_GU(solid,P3P)] 

    ex = FreeCAD.Vector(1,0,0)
    ey = FreeCAD.Vector(0,1,0)
    ez = FreeCAD.Vector(0,0,1)    
    Surfaces=UF.Surfaces_dict()    

    for solid_GU in solidParts:
       flag_inv = CD.isInverted(solid_GU.solid)
     
       # Get the parameter of all faces of the solid
       # Add auxillary planes for Cylinders and Cones                    
       for face in solid_GU.Faces:
         surf = str(face.Surface)
         if surf == '<Plane object>' and kind in ['Planes','All'] :
            normal = face.Surface.Axis
            pos=face.CenterOfMass
            dim1 = face.ParameterRange[1]-face.ParameterRange[0]
            dim2 = face.ParameterRange[3]-face.ParameterRange[2]
            plane = UF.GEOUNED_Surface(('Plane',(pos,normal,dim1,dim2)),UniverseBox)
            if MakeObj : plane.buildSurface()
            Surfaces.addPlane(plane,fuzzy)
    
         elif surf == '<Cylinder object>' :
           dir  = face.Surface.Axis
           orig = face.Surface.Center
           rad  = face.Surface.Radius
           dimL = face.ParameterRange[3]-face.ParameterRange[2]
           if kind in ['Cyl','All'] :
             cylinder = UF.GEOUNED_Surface(('Cylinder',(orig,dir,rad,dimL)),UniverseBox)
             if MakeObj : cylinder.buildSurface()
             Surfaces.addCylinder(cylinder,fuzzy)

           if kind in ['Planes','All'] :
             #add planes if cylinder axis is cut by a plane (plane quasi perpedicular to axis)
             for p in CylBoundPlanes(face,UniverseBox) :
                 if MakeObj : p.buildSurface()
                 Surfaces.addPlane(p)

         elif ( surf == '<Cone object>') :
           dir = face.Surface.Axis
           apex = face.Surface.Apex
           halfAngle = face.Surface.SemiAngle         
           dimL = face.ParameterRange[3]-face.ParameterRange[2]
           dimR = face.Surface.Radius
           if kind in ['Cone','All'] :
             cone = UF.GEOUNED_Surface(('Cone',(apex,dir,halfAngle,dimL,dimR)),UniverseBox)
             if MakeObj : cone.buildSurface()
             Surfaces.addCone(cone)
        
           if kind in ['Planes','All'] :
             for p in CylBoundPlanes(face,UniverseBox) :
                 if MakeObj : p.buildSurface()
                 Surfaces.addPlane(p)

         elif surf[0:6] == 'Sphere' and kind in ['Sph','All']:
           rad = face.Surface.Radius
           pnt = face.Surface.Center
           sphere = UF.GEOUNED_Surface(('Sphere',(pnt,rad)),UniverseBox)
           if MakeObj : sphere.buildSurface() 
           Surfaces.addSphere(sphere)
         
           if kind in ['Planes','All'] :
             for p in CylBoundPlanes(face,UniverseBox) :
                 if MakeObj : p.buildSurface()
                 Surfaces.addPlane(p)

         elif (surf == '<Toroid object>'):
           if  kind in ['Tor','All']:
              radMaj = face.Surface.MajorRadius
              radMin = face.Surface.MinorRadius
              center = face.Surface.Center
              dir    = face.Surface.Axis
              torus = UF.GEOUNED_Surface(('Torus',(center,dir,radMaj,radMin)),UniverseBox)
              if MakeObj : torus.buildSurface() 
              Surfaces.addTorus(torus)

           if kind in ['Planes','All'] :
              for p in TorusBoundPlanes(face,UniverseBox) :
                 if MakeObj : p.buildSurface()
                 Surfaces.addPlane(p)

         elif surf == '<Plane object>' and kind == 'Plane3Pts' :
            pos=face.CenterOfMass
            normal = face.Surface.Axis
            dim1 = face.ParameterRange[1]-face.ParameterRange[0]
            dim2 = face.ParameterRange[3]-face.ParameterRange[2]
            points = tuple(v.Point for v in face.Vertexes)

            plane = UF.GEOUNED_Surface(('Plane3Pts',(pos,normal,dim1,dim2,points)),UniverseBox)
            if MakeObj : plane.buildSurface()
            Surfaces.addPlane(plane,fuzzy)
         
    return Surfaces                  
 
    
def isAlreadyInPlanes(plane,Array):

    for elem in Array:
        if (plane.Surface.Axis.cross(elem.Surface.Axis)==FreeCAD.Vector(0,0,0) and isInPlane(plane.Surface.Position,elem.Surface)):
            return True
    
    return False
    
#
#
#   Check if to faces are joint
#    
def ContiguousFace(face1,face2):
    return face1.distToShape(face2)[0] < tol.distance
    
def SameFaces(Faces):
    Connection=OrderedDict()
    if len(Faces)==1:
        return []
    
    for i, face1 in enumerate(Faces):
        Couples=[]
        if not Faces[i+1:]:
            continue
        for j, face2 in enumerate(Faces[i+1:]):
            if ContiguousFace(face1,face2):
                
                Couples.append(i+1+j)
                
        Connection[i]=Couples
                      
    lista=Connection[0]
    Connection.popitem(0)
    
    if len(Connection)==0: # solo los elementos de la lista de la superficie 0
        return lista
    
    if not lista: # ninguna face está conecta conecta con la superficie 0
        return lista
    
    for elem in Connection:
        if elem in lista: # que la key esta en lista implica que sus dependencias estan
            lista.extend(Connection[elem])
        else:
            for elem2 in Connection[elem]:
                if elem2 in lista: # si una de sus dependencias esta en lista lo esta la clave
                    lista.append(elem)

    return list(set(lista))
    
# Tolerance in this function are not the general once
# function should be reviewed
def GenPlaneCylinder(face,solid):
    Surf=face.Surface
    rad=Surf.Radius
    if face.Area < tol.min_area_gen:
        return None
    
    UVNodes=[] 
    face_index_0=[solid.Faces.index(face)]

    try:
        face.tessellate(0.1)
        UVNodes.append(face.getUVNodes())
    except RuntimeError:
        PR=face.ParameterRange
        UVNode1=(PR[0],PR[2])
        UVNode2=(PR[1],PR[3])
        UVNodes.append([UVNode1,UVNode2])
        
    
    for i, face2 in enumerate(solid.Faces):
              
        if str(face2.Surface) == '<Cylinder object>' and not(face2.isEqual(face)):
            if (face2.Surface.Axis.isEqual(face.Surface.Axis,1e-5) and face2.Surface.Radius == rad and isInLine(face2.Surface.Center,face.Surface.Axis,face.Surface.Center)):
                face_index_0.append(i)
                            
    # prueba SameFaces, parece ok
    Faces_p=[]
    for ind in face_index_0:
        Faces_p.append(solid.Faces[ind])
    
    face_index=[face_index_0[0]] # la face de entrada

    for k in SameFaces(Faces_p):
        face_index.append(face_index_0[k])
    
    # commented to let plane cut rounded corner 
    # if len(face_index_0)==len(face_index): #Esto evita un corte creo!!! no debería ser así en conversion
    #    return None
    
    for ind in reversed(face_index[1:]):
        if solid.Faces[ind].Area <= tol.min_area_gen :
            face_index.remove(ind)

        else:
            face2 = solid.Faces[ind]
            try:
                face2.tessellate(0.1)
                UVNodes.append(face2.getUVNodes())
            except RuntimeError:
                PR=face2.ParameterRange
                UVNode1=(PR[0],PR[2])
                UVNode2=(PR[1],PR[3])
                UVNodes.append([UVNode1,UVNode2])
    
    AngleRange=0.0
    
    Uval=[]
    for index in face_index:
        Range=solid.Faces[index].ParameterRange
        AngleRange=AngleRange+abs(Range[1]-Range[0])
        #if not(Range[0] in Uval) and not(Range[1] in Uval): 
        Uval.append(Range[0])
        Uval.append(Range[1])
             
    if (twoPi-AngleRange < 1.e-2 or AngleRange < 1.e-2):
        return None
    
    Uval_str_cl=[]
    
    for i,elem1 in enumerate(Uval):
        num_str1='%11.4E' %elem1
        if abs(elem1)<1.0e-5:
            num_str1='%11.4E' %0.0
    
        if not(isDuplicateInList(num_str1,i,Uval)):
            Uval_str_cl.append(num_str1)
    
    if len(Uval_str_cl) < 2:
        if opt.verbose : print ('GenPlaneCylinder : Uval_str_cl should no be void ')
        return None

    face_index_2=[face_index[0],face_index[0]]
    
    Node_min=UVNodes[0][0]
    Node_max=UVNodes[0][1]
  
    dif1_0=abs(float(Uval_str_cl[0])-Node_min[0])
    dif2_0=abs(float(Uval_str_cl[1])-Node_max[0])
           
    # searching for minimum and maximum angle points
   
    for j,Nodes in enumerate(UVNodes): 
        for elem in Nodes:
            val='%11.4E' %elem[0]
            dif1=abs(float(Uval_str_cl[0])-elem[0])
            dif2=abs(float(Uval_str_cl[1])-elem[0])

            if abs(elem[0])<1.0e-5:
                val='%11.4E' %0.0
            if (dif1<dif1_0):
                Node_min=elem
                face_index_2[0]=face_index[j]
                dif1_0=dif1
            if (dif2<dif2_0):
                Node_max=elem
                face_index_2[1]=face_index[j]
                dif2_0=dif2           
    
    V1=solid.Faces[face_index_2[0]].valueAt(Node_min[0],Node_min[1])
    V2=solid.Faces[face_index_2[1]].valueAt(Node_max[0],Node_max[1])
    
    if (V1.isEqual(V2,1e-5)):
        if opt.verbose : print('Error in the additional plane definition')
        return None
    
    normal=V2.sub(V1).cross(face.Surface.Axis)
    
    plane=Part.Plane(V1,normal).toShape()
    
    return plane

# Tolerance in this function are not the general once
# function should be reviewed    
def GenPlaneCone(face,solid):

    Surf=face.Surface
    rad=Surf.Radius
    
    if face.Area < tol.min_area_gen :
        return None
           
    UVNodes=[] 
    face_index_0=[solid.Faces.index(face)]

    try:
        face.tessellate(0.1)
        UVNodes.append(face.getUVNodes())
    except RuntimeError:
        PR=face.ParameterRange
        UVNode1=(PR[0],PR[2])
        UVNode2=(PR[1],PR[3])
        UVNodes.append([UVNode1,UVNode2])

    for i, face2 in enumerate(solid.Faces):
              
        if str(face2.Surface) == '<Cone object>' and not(face2.isEqual(face)):
                        
            if (face2.Surface.Axis.isEqual(face.Surface.Axis,1e-5) and face2.Surface.Apex.isEqual(face.Surface.Apex,1e-5) and (face2.Surface.SemiAngle-face.Surface.SemiAngle) < 1e-6):

                face_index_0.append(i)
    
    Faces_p=[]
    for ind in face_index_0:
        Faces_p.append(solid.Faces[ind])
    
    face_index=[face_index_0[0]] # la face de entrada
    
    for k in SameFaces(Faces_p):
        face_index.append(face_index_0[k])
     
    # same as cylinder commennt     
    #if len(face_index_0)==len(face_index):
    #    return None
    
    for ind in reversed(face_index[1:]):
        if solid.Faces[ind].Area<=tol.min_area_gen:
            face_index.remove(ind)
        else:
            face2 = solid.Faces[ind]
            try:
                face2.tessellate(0.1)
                UVNodes.append(face2.getUVNodes())
            except RuntimeError:
                PR=face2.ParameterRange
                UVNode1=(PR[0],PR[2])
                UVNode2=(PR[1],PR[3])
                UVNodes.append([UVNode1,UVNode2])
    
    AngleRange=0.0
    
    Uval=[]
     
    for index in face_index:
        Range=solid.Faces[index].ParameterRange
        AngleRange=AngleRange+abs(Range[1]-Range[0])
        Uval.append(Range[0])
        Uval.append(Range[1])
    if (twoPi-AngleRange < 1.e-2 or AngleRange<1e-2):
        return None
    
    Uval_str_cl=[]
    
    for i,elem1 in enumerate(Uval):
        num_str1='%11.4E' %elem1
        if abs(elem1)<1.0e-5:
            num_str1='%11.4E' %0.0
        if not(isDuplicateInList(num_str1,i,Uval)):
            Uval_str_cl.append(num_str1)
    
    if len(Uval_str_cl) < 2:
        if opt.verbose : print ('GenPlaneCone : Uval_str_cl should no be void ')
        return None

    face_index_2=[face_index[0],face_index[0]]
    
    Node_min=UVNodes[0][0]
    Node_max=UVNodes[0][1]
    dif1_0=abs(float(Uval_str_cl[0])-Node_min[0])
    dif2_0=abs(float(Uval_str_cl[1])-Node_max[0])
           
    # searching for minimum and maximum angle points
    for j,Nodes in enumerate(UVNodes): 
        for elem in Nodes:
            val='%11.4E' %elem[0]
            dif1=abs(float(Uval_str_cl[0])-elem[0])
            dif2=abs(float(Uval_str_cl[1])-elem[0])

            if abs(elem[0])<1.0e-5:
                val='%11.4E' %0.0
            if (dif1<dif1_0):
                Node_min=elem
                face_index_2[0]=face_index[j]
                dif1_0=dif1
            if (dif2<dif2_0):
                Node_max=elem
                face_index_2[1]=face_index[j]
                dif2_0=dif2
    
    V1=solid.Faces[face_index_2[0]].valueAt(Node_min[0],Node_min[1])
    V2=solid.Faces[face_index_2[1]].valueAt(Node_max[0],Node_max[1])
    
    if (V1.isEqual(V2,1e-5)):
        if opt.verbose : print('Error in the additional plane definition')
        return None
    
    #normal=V2.sub(V1).cross(face.Surface.Axis)
    
    #plane=Part.Plane(V1,normal).toShape()
    plane=Part.Plane(V1,V2,face.Surface.Apex).toShape()
    
    return plane
    
def Plane2ndOrder(solid_GU,face,flag_inv,convex=True):
    planes = []

    if face is None :
       for face in solid_GU.Faces:
          if (flag_inv):
              orient_temp=face.Orientation
              if orient_temp=='Forward':
                  orient='Reversed'
              elif orient_temp=='Reversed':
                  orient='Forward'
          else:
              orient=face.Orientation
    
          if convex and orient=='Reversed': # convex face condition
              pp = None
              if str(face.Surface) == '<Cylinder object>': 
                  pp = GenPlaneCylinder(face,solid_GU)
              elif str(face.Surface) == '<Cone object>':
                  pp = GenPlaneCone(face,solid_GU)
              if pp is not None : 
                planes.append(pp)
            #elif str(face.Surface)[0:6] == 'Sphere':
            #    return GenPlaneSphere(face,solid_GU) 

    else:
       if (flag_inv):
           orient_temp=face.Orientation
           if orient_temp=='Forward':
               orient='Reversed'
           elif orient_temp=='Reversed':
               orient='Forward'
       else:
           orient=face.Orientation

       if not convex and orient=='Forward':
            if str(face.Surface) == '<Cylinder object>': 
               pp = GenPlaneCylinder(face,solid_GU)
            elif str(face.Surface) == '<Cone object>':
               pp = GenPlaneCone(face,solid_GU)
            if pp is not None : 
               planes.append(pp)

    return planes

def SplitPlanes(Solids,UniverseBox,newVersion=True):
    if newVersion :
        return SplitPlanes_new(Solids,UniverseBox)
    else:
        return SplitPlanes_org(Solids,UniverseBox)

def SplitPlanes_new(Solids,UniverseBox):
    Bases = Solids[:]
    simpleSolid = []

    # Cut with real planes defining solid faces
    while True: 
      newBases = []
      for base in Bases:
         cutSolids = splitPPlanes_new(base,UniverseBox)
         if len(cutSolids) == 1:
           simpleSolid.extend(cutSolids)
         else:
           newBases.extend(cutSolids)
      if len(newBases) == 0 :
         break
      else:
         Bases = newBases

    return simpleSolid,0

 
def SplitPlanes_org(Solids,UniverseBox):
    Bases = []
    icount=0
    err =0
    for sol in Solids :
      Bases.append((sol,[]))
       
# Loop  1 2 3   1 2 3  1 2 3   1 2 3   1 2 3  1 2 3
# order x y z   x z y  y x z   y z x   z x y  z y x
# index 0 0 0   0 1 0  1 0 0   1 1 0   2 0 0  2 1 0

    # planes orthogonal to X,Y,Z axis
    for iplane in range(3) :
       newBase =[]
       for item in Bases :
           
         base  = item[0]
         index = item[1]
         SPlanes = ExtractSurfaces(base,'Planes',UniverseBox,MakeObj=True)
         Planes = [SPlanes['PX'],SPlanes['PY'],SPlanes['PZ']]
         for i in index :
            del(Planes[i])
         pmin = len(Planes[0])
         imin = 0                
         for i in range(2-iplane): 
           if pmin > len(Planes[i+1]) :
             pmin = len(Planes[i+1])
             imin = i+1  
         
         if len(Planes[imin]) != 0:
            Tools = []
            for p in Planes[imin] :
              p.buildSurface() 
              Tools.append(p.shape)
            comsolid=BOPTools.SplitAPI.slice(base,Tools,'Split', tolerance = opt.splitTolerance)
            if len(comsolid.Solids) == 1 :
               if abs(comsolid.Solids[0].Volume-base.Volume)/base.Volume > tol.relativePrecision :
                  if opt.verbose : print('Warning. Part of the split object is missing original base is used instead',\
                        abs(comsolid.Solids[0].Volume-base.Volume)/base.Volume,comsolid.Solids[0].Volume,base.Volume)
                  base.exportStep('tmp_base.stp')
                  s = Part.Shape()
                  s.read('tmp_base.stp')                  
                  solid_list = s.Solids
                  err=1
               else:
                  solid_list = comsolid.Solids
            elif len(comsolid.Solids)> 1 :
               solid_list = comsolid.Solids
            else:
               solid_list = [base]  
         else :
            solid_list = [base]

         newindex = index[:]
         newindex.append(imin)
         for sol in solid_list:
           newBase.append((sol,newindex))
       Bases = newBase

    XYZBases = []
    for ii,s in enumerate(Bases):
      XYZBases.append(s[0])

    simpleSolid = []
    Bases = XYZBases

    # other planes
    while True: 
      newBases = []
      for base in Bases:
         cutSolids = splitPPlanes_org(base,UniverseBox)
         if len(cutSolids) == 1:
           simpleSolid.extend(cutSolids)
         else:
           newBases.extend(cutSolids)
      if len(newBases) == 0 :
         break
      else:
         Bases = newBases
    return simpleSolid,err


class parallelPlanes:
    def __init__(self,plane):
       self.Axis = plane.Surf.Axis
       self.elements = [plane]
       self.count = 1

    def append(self,plane):
       if isParallel(plane.Surf.Axis,self.Axis,0.1):
          self.elements.append(plane)
          self.count += 1
          return True
       else: 
          return False 

    def sort(self):
       pos=[]
       for i,p in enumerate(self.elements):
          d = self.Axis.dot(p.Surf.Position)  
          pos.append((d,i))
       pos.sort()

       sort = []
       for p in pos:
         sort.append(self.elements[p[1]])
       self.elements = sort

   
def sortPlanes(PlaneList,sortElements=False):
    if not PlaneList : return []
    newList = [parallelPlanes(PlaneList[0])]
    for p in PlaneList[1:]:
       found = False
       for pp in newList:
         if pp.append(p) : 
            found = True
            break
       if not found :
          newList.append(parallelPlanes(p)) 

    lenList = []
    for i,pp in enumerate(newList):
      if sortElements : pp.sort()
      lenList.append((pp.count,i))
    lenList.sort()

    sortedPlanes=[]
    for l in lenList:
       sortedPlanes.append(newList[l[1]].elements)
          
    return sortedPlanes
         
def splitPPlanes_new(solid,UniverseBox):
    SPlanes = ExtractSurfaces(solid,'Planes',UniverseBox)
  
    Planes = [] 
    for P in ('PX','PY','PZ','P'):
       Planes.extend(SPlanes[P])

    Planes = sortPlanes(Planes)
    
    if len(Planes) == 0 : return [solid]
    if len(Planes[-1]) < opt.nPlaneReverse : Planes.reverse()
    outSolid = [solid]
    for pp in Planes :
       for p in pp : p.buildSurface()
       tools = tuple( p.shape for p in pp) 
       comsolid=BOPTools.SplitAPI.slice(solid,tools,'Split',tolerance =opt.splitTolerance)         

       if len(comsolid.Solids) > 1 :
         outSolid = comsolid.Solids
         break
    return outSolid   
  
def splitPPlanes_org(solid,UniverseBox):
    SPlanes = ExtractSurfaces(solid,'Planes',UniverseBox)
    
    if len(SPlanes['P']) == 0 : return [solid]
    outSolid = [solid]
    for p in SPlanes['P'] :
       p.buildSurface()
       comsolid=BOPTools.SplitAPI.slice(solid,[p.shape],'Split', tolerance = opt.splitTolerance)
       if len(comsolid.Solids) > 1 :
         outSolid = comsolid.Solids
         break
    return outSolid   
        
def Split2ndOrder(Solids,UniverseBox):
    err = 0
    Base = Solids
    for kind in ['Cyl','Cone','Sph','Tor'] :                    
      kindBase = []
      while True :    
        cutBase=[]
        for solid in Base :
           Surfaces = ExtractSurfaces(solid,kind,UniverseBox)
           if len(Surfaces[kind]) == 0 :
              kindBase.append(solid)
           else:    
              simple = True
              for s in Surfaces[kind]:
                 s.buildSurface()
                 try :
                    comsolid=BOPTools.SplitAPI.slice(solid,[s.shape],'Split', tolerance = opt.splitTolerance)

                    solidsInCom = []
                    for s in comsolid.Solids :
                        if s.Volume > tol.min_volume : solidsInCom.append(s)
                        else:
                            if opt.verbose : 
                                print('Warning: splitComponent degenerated solids are produced', \
                                        part.Volume)

                        
                    if len(solidsInCom) > 1 :
                       simple =False
                       cutBase.extend(solidsInCom)
                       break 
                 except :
                    if opt.verbose : print('Failed split base with {} surface'.format(s.shape.Faces[0].Surface))
                    err += 2

              if simple :         
                 kindBase.append(solid)

        if len(cutBase) == 0 :
           Base = kindBase
           break
        else:
           Base = cutBase

    return Base,err
    
def Split2ndOrderPlanes(Solids):
    err = 0
    simpleSolid = []
    Bases = Solids
    while True:
        newBases = []
        for base in Bases:
            cutSolids,err = split2ndOPlane(base)
            if len(cutSolids) == 1:
                simpleSolid.extend(cutSolids)
            else:
                newBases.extend(cutSolids)
        if len(newBases) == 0 :
            break
        else:
            Bases = newBases
         
    return simpleSolid,err
    
def split2ndOPlane(solid):

    Planes=[]
    err = 0
    flag_inv=CD.isInverted(solid)
    solid_GU=GU.solid_GU(solid)
    planes = Plane2ndOrder(solid_GU,None,flag_inv,convex=True)
    if not planes  : return [solid],err

    for p in planes:
      comsolid=BOPTools.SplitAPI.slice(solid,[p],'Split', tolerance = opt.splitTolerance)
      if not comsolid.Solids:
         comsolid=solid
         continue
      if len(comsolid.Solids) > 1 :
        removed, err = removeSolids(comsolid.Solids)
        comsolid=Part.makeCompound(removed)
        if len(removed) > 1: break

    return comsolid.Solids,err
    
def removeSolids(Solids):
    err = 0
    Solids_Clean=[]
    for solid in Solids:
        if solid.Volume <= tol.min_volume:
            if opt.verbose : print('Warning: removeSolids degenerated solids are produced',solid.Volume)
            err=2
            continue
        Solids_Clean.append(solid)
    
    return Solids_Clean,err
                        
def splitComponent(solidShape,UniverseBox):
    err  = 0
    err2 = 0

    Volini=solidShape.Volume
    Solids=solidShape.Solids
    # Split with explicit planes bounding the solid and
    # implicit planes interface of two 2nd order surfaces
    split0,err = SplitPlanes(Solids,UniverseBox,opt.newSplitPlane)
    # Split with explicit 2nd order surfaces bounding the solid
    
    split1,err1 = Split2ndOrder(split0,UniverseBox)
    err += err1

    split,err2  = Split2ndOrderPlanes(split1)
    err += err2
    Pieces = []
    for part in split:
       if abs(part.Volume) <= tol.min_volume:
          if opt.verbose : print('Warning: splitComponent degenerated solids are produced',part.Volume)
          err += 2
          continue
       Pieces.append(part)                  
       
    comsolid=Part.makeCompound(Pieces)
 
    Volend=comsolid.Volume
    Volch=(Volini-Volend)/Volini
 
    if (abs(Volch)>1e-2): # 1% of Volume change
        if opt.verbose : print('Warning: Volume has changed after decomposition: %11.4E' %Volch)
        err += 4

    return comsolid, err
    
def SplitSolid(solidShape,UniverseBox):

    solidParts = []

    for solid in solidShape.Solids :
       
       explode = splitFullCylinder(solid)
       piece, err = splitComponent(explode,UniverseBox)
       solidParts.append(piece)
          
    
    return Part.makeCompound(solidParts) ,err   
    
