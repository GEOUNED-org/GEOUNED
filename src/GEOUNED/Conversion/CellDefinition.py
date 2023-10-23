############################
# Module for Cell definiton #
#############################
import FreeCAD, Part
from   GEOUNED.Utils.Functions import GEOUNED_Surface
from   GEOUNED.Utils.BasicFunctions_part1  \
       import isParallel, isOposite, isInLine, signPlane, isSameValue
from    GEOUNED.Utils.booleanFunction import BoolSequence

import  GEOUNED.Utils.BasicFunctions_part2 as BF 
import GEOUNED.Utils.Qform as QF
import GEOUNED.Utils.Geometry_GU as GU
from GEOUNED.Utils.Options.Classes import Tolerances as tol
from GEOUNED.Utils.Options.Classes import Options as opt
import math
from collections import OrderedDict

def getId(facein, Surfaces):
    
    surfin = str(facein)
    if surfin == '<Plane object>' :
        if   isParallel(facein.Axis,FreeCAD.Vector(1,0,0),tol.pln_angle) :
          P = 'PX'
        elif isParallel(facein.Axis,FreeCAD.Vector(0,1,0),tol.pln_angle) :
          P = 'PY'
        elif isParallel(facein.Axis,FreeCAD.Vector(0,0,1),tol.pln_angle) :
          P = 'PZ'
        else:
          P = 'P'

        for s in Surfaces[P] :
          if BF.isSamePlane(facein,s.Surf,dtol=tol.pln_distance,atol=tol.pln_angle,relTol=tol.relativeTol) :
              return s.Index
    
    elif surfin == '<Cylinder object>' :
        for s in Surfaces['Cyl'] :
           if BF.isSameCylinder(facein,s.Surf,dtol=tol.cyl_distance,atol=tol.cyl_angle,relTol=tol.relativeTol) :  
              return s.Index 
        
    elif surfin == '<Cone object>' :
        for s in Surfaces['Cone'] :
           if BF.isSameCone(facein,s.Surf,dtol=tol.kne_distance,atol=tol.kne_angle,relTol=tol.relativeTol):
              return s.Index
         
    elif surfin[0:6] == 'Sphere':
        for s in Surfaces['Sph'] :
           if BF.isSameSphere(facein,s.Surf,tol.sph_distance,relTol=tol.relativeTol):
              return s.Index

    elif surfin == '<Toroid object>' :
        for s in Surfaces['Tor'] :
           if BF.isSameTorus(facein,s.Surf,dtol=tol.tor_distance,atol=tol.tor_angle,relTol=tol.relativeTol):
              return s.Index
                
    return 0
   
def isInverted(solid):

    face=solid.Faces[0]

    #u=(face.Surface.bounds()[0]+face.Surface.bounds()[1])/2.0 # entre 0 y 2pi si es completo
    #v=face.Surface.bounds()[0]+(face.Surface.bounds()[3]-face.Surface.bounds()[2])/3.0 # a lo largo del eje
    Range=face.ParameterRange
    u=(Range[1]+Range[0])/2.0
    v=(Range[3]+Range[2])/2.0
    
    if (str(face.Surface)=='<Cylinder object>'):
        dist1=face.Surface.value(u,v).distanceToLine(face.Surface.Center,face.Surface.Axis)
        dist2=face.Surface.value(u,v).add(face.Surface.normal(u,v).multiply(1.0e-6)).distanceToLine(face.Surface.Center,face.Surface.Axis)
        if ((dist2-dist1)<0.0):
            # The normal of the cylinder is going inside
            return True
    elif (str(face.Surface)=='<Cone object>' ):
        dist1=face.Surface.value(u,v).distanceToLine(face.Surface.Apex,face.Surface.Axis)
        dist2=face.Surface.value(u,v).add(face.Surface.normal(u,v).multiply(1.0e-6)).distanceToLine(face.Surface.Apex,face.Surface.Axis)
        if ((dist2-dist1)<0.0):
            # The normal of the cylinder is going inside
            return True
    # MIO
    elif str(face.Surface)[0:6]=='Sphere':
        # radii = point - center
        radii  = face.Surface.value(u,v).add(face.Surface.Center.multiply(-1))
        radiiB  = face.Surface.value(u,v).add(face.Surface.normal(u,v).multiply(1.0e-6)).add(face.Surface.Center.multiply(-1))
        #radiiB  = radii.add( face.Surface.normal(u,v).multiply(1.0e-6) )
        if ((radiiB.Length-radii.Length)<0.0):
            # An increasing of the radii vector in the normal direction decreases the radii: oposite normal direction
            return True
           
    elif  (str(face.Surface)=='<Plane object>'):
        dist1=face.CenterOfMass.distanceToPoint(solid.BoundBox.Center)
        dist2=face.CenterOfMass.add(face.normalAt(u,v).multiply(1.0e-6)).distanceToPoint(solid.BoundBox.Center)
        point2=face.CenterOfMass.add(face.normalAt(u,v).multiply(1.0e-6))
        if (solid.isInside(point2,1e-7,False)):
            return True
            
    return False
    

def GenPlane(face,solid):

    """ Generate an additional plane when convex surfaces of second order are presented as a face of the solid """

    
    surf = face.Surface
    if (str(surf) == '<Cylinder object>'):
        return GenPlaneCylinder(face,solid)
    if (str(surf) == '<Cone object>'):
        return GenPlaneCone(face,solid)
    
    rad=surf.Radius
    axis=surf.Axis
    pnt=surf.Center

    Same_Faces=[]
    Same_Faces.append(face)
    
    for face2 in solid.Faces:
        
# modif Pat 03/07/2023 
# here face cannot be either cone nor cylinder so if face2 is cone or cylinder it cannot be equal to face

#        if str(face2.Surface) == '<Cone object>' and not(face2.isEqual(face)):
#            if (face2.Surface.Axis.isEqual(face.Surface.Axis,1e-5) and face2.Surface.Apex.isEqual(face.Surface.Center,1e-5) and (face2.Surface.SemiAngle-face.Surface.SemiAngle)<1e-2):
                #print 'Warning: coincident cone faces are the same'
#                Same_Faces.append(face2)
        
#        if str(face2.Surface) == '<Cylinder object>' and not(face2.isEqual(face)):
#            if (face2.Surface.Axis.isEqual(face.Surface.Axis,1e-5) and face2.Surface.Radius == rad and isInLine(face2.Surface.Center,face.Surface.Axis,face.Surface.Center)):
                #print 'Warning: coincident cylinder faces are the same'
#                Same_Faces.append(face2)
                
        if str(surf)[0:6] == 'Sphere' and not(face2.isEqual(face)):
            if (face2.Surface.Center == face.Surface.Center and face2.Surface.Radius == face.Surface.Radius):
                #print 'Warning: coincident sphere faces are the same'
                Same_Faces.append(face2)

    Vertexes=[]
    #print Same_Faces
    for face3 in Same_Faces:
        for Vertex in face3.Vertexes:
            Vertexes.append(Vertex.Point)
            
    Vectors=[] # free vectors, right directed but not located
    Vectors2=[] # right coordinates from origin 
    if (str(surf) == '<Cylinder object>' or str(surf) == '<Cone object>'):
        for Point in Vertexes:
            Vec=pnt.sub(Point).projectToPlane(pnt,axis)
            Vec2=Point.projectToPlane(pnt,axis)
            Vectors.append(Vec) # free vectors
            Vectors2.append(Vec2)
    elif (str(surf)[0:6] == 'Sphere'):
        for Point in Vertexes:
            Vec=pnt.sub(Point)
            Vectors.append(Vec)
            
    angle_max=0.0
    for i,V1 in enumerate(Vectors):
        for j,V2 in enumerate(Vectors):
            angle=V1.getAngle(V2)
            if (angle > angle_max):
                Va=V1
                Vb=V2
                Va2=Vectors2[i]
                Vb2=Vectors2[j]
                angle_max=angle
                
    if (str(surf) == '<Cylinder object>'):
        normal=Va.sub(Vb).cross(surf.Axis)
    elif (str(surf) == '<Cone object>'):
        normal=Va.sub(Vb).cross(Va2.sub(surf.Apex))
    
    plane=Part.makePlane(rad,rad,pnt,normal)
    vec_on_plane = plane.Vertexes[3].Point.sub(plane.Vertexes[0].Point)
    new_pos = plane.Vertexes[0].Point.sub(vec_on_plane.multiply(10.0))
    plane_center = Part.makePlane(rad*100,rad*100,new_pos,normal)
    plane=plane_center
    
    sign= signPlane(face.Vertexes[0].Point,plane)
    normal.multiply(float(sign))
    
    plane=Part.Plane(pnt,normal).toShape()
    
    if (solid.distToShape(plane)[0]<1e-8):
        return None
    else:
        return Part.Plane(Va2,normal).toShape()
    
def GenPlaneCylinder(face,solid):

    Surf=face.Surface
    rad=Surf.Radius

    if (str(Surf)!='<Cylinder object>'):
        return None
        
    UVNodes=[]
    face_index=[solid.Faces.index(face)]
    
    try:    
        UVNodes.append(face.getUVNodes())
    except RuntimeError:
        face.tessellate(100.0)
        UVNodes.append(face.getUVNodes())
        

    for i, face2 in enumerate(solid.Faces):
              
        if face2.Area < 1e-2 : continue 
        if str(face2.Surface) == '<Cylinder object>' and not(face2.isEqual(face)):
            if (face2.Surface.Axis.isEqual(face.Surface.Axis,1e-5) and face2.Surface.Radius == rad and isInLine(face2.Surface.Center,face.Surface.Axis,face.Surface.Center)):
                #print 'Warning: coincident cylinder faces are the same'
                try:    
                    UVNodes.append(face2.getUVNodes())
                except RuntimeError:
                    face2.tessellate(100.0)
                    UVNodes.append(face2.getUVNodes())
                face_index.append(i)
    
    AngleRange=0.0
    
    Uval=[]
    
    for index in face_index:
        Range=solid.Faces[index].ParameterRange
        AngleRange=AngleRange+abs(Range[1]-Range[0])
        if not(Range[0] in Uval) and not(Range[1] in Uval): 
            Uval.append(Range[0])
            Uval.append(Range[1])
            
    if (2.0*math.pi-AngleRange<1e-2):
        return None
    
    Uval_str_cl=[]
    for i,elem1 in enumerate(Uval):
        num_str1='%11.4E' %elem1
        if abs(elem1)<1.0e-5:
            num_str1='%11.4E' %0.0
        if not(BF.isDuplicateInList(num_str1,i,Uval)):
            Uval_str_cl.append(num_str1)
    
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
        if opt.verbose : print('Error in the additional place definition')
        return None
        
    normal=V2.sub(V1).cross(face.Surface.Axis)
    plane=Part.Plane(V1,normal).toShape()
    
    return plane
    
def GenPlaneCone(face,solid):

    Surf=face.Surface
    rad=Surf.Radius
    Axis=face.Surface.Axis
    if (str(Surf)!='<Cone object>'):
        return None
        
    UVNodes=[]
    face_index=[]

    try:    
        UVNodes.append(face.getUVNodes())
    except RuntimeError:
        face.tessellate(100.0)
        UVNodes.append(face.getUVNodes())
        
    face_index.append(solid.Faces.index(face))

    for i, face2 in enumerate(solid.Faces):
        if face2.Area < 1e-2 : continue 
              
        if str(face2.Surface) == '<Cone object>' and not(face2.isEqual(face)):
                        
            if (face2.Surface.Axis.isEqual(face.Surface.Axis,1e-5) and face2.Surface.Apex.isEqual(face.Surface.Apex,1e-5) and (face2.Surface.SemiAngle-face.Surface.SemiAngle) < 1e-6):
                try:    
                    UVNodes.append(face2.getUVNodes())
                except RuntimeError:
                    face2.tessellate(100.0)
                    UVNodes.append(face2.getUVNodes())
                face_index.append(i)
    
    AngleRange=0.0
    
    Uval=[]
    
    for index in face_index:
        Range=solid.Faces[index].ParameterRange
        AngleRange=AngleRange+abs(Range[1]-Range[0])
        Uval.append(Range[0])
        Uval.append(Range[1])
                
    if (2.0*math.pi-AngleRange<1e-2):
        return None
    
    Uval_str_cl=[]
    
    for i,elem1 in enumerate(Uval):
        num_str1='%11.4E' %elem1
        if abs(elem1)<1.0e-5:
            num_str1='%11.4E' %0.0
        if not(BF.isDuplicateInList(num_str1,i,Uval)):
            Uval_str_cl.append(num_str1)
    
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
        if opt.verbose : print('Error in the additional place definition')
        return None
       
    plane=Part.Plane(V1,V2,face.Surface.Apex).toShape()
    
    return plane


def GenTorusAnnexUPlanes(face,Uparams):
    if   isParallel(face.Surface.Axis,FreeCAD.Vector(1,0,0),tol.tor_angle) :
        axis = FreeCAD.Vector(1,0,0)
    elif isParallel(face.Surface.Axis,FreeCAD.Vector(0,1,0),tol.tor_angle) :
        axis = FreeCAD.Vector(0,1,0)
    elif isParallel(face.Surface.Axis,FreeCAD.Vector(0,0,1),tol.tor_angle) :
        axis = FreeCAD.Vector(0,0,1)

    center = face.Surface.Center
    p1 = face.valueAt(Uparams[0],0.)
    p2 = face.valueAt(Uparams[1],0.)
    pmid = face.valueAt(0.5*(Uparams[0]+Uparams[1]),0.)
        
    if isSameValue(abs(Uparams[1]-Uparams[0]),math.pi,tol.value):
        d = axis.cross(p2-p1)
        d.normalize()
        if pmid.dot(d) < 0:  d = -d
        return ((center,d,face.Surface.MajorRadius),tuple()),False

    else :
        d1 = axis.cross(p1)
        d1.normalize()
        if pmid.dot(d1) < 0:  d1 = -d1

        d2 = axis.cross(p2)
        d2.normalize()
        if pmid.dot(d2) < 0:  d2 = -d2

        if Uparams[1]-Uparams[0] < math.pi :
          return ((center,d1,face.Surface.MajorRadius,face.Surface.MajorRadius),(center,d2,face.Surface.MajorRadius,face.Surface.MajorRadius)),False # ( d1 d2 )
        else:
          return ((center,d1,face.Surface.MajorRadius,face.Surface.MajorRadius),(center,d2,face.Surface.MajorRadius,face.Surface.MajorRadius)),True  # (d1 : d2)
        
def GenTorusAnnexVSurface(face,Vparams,forceCylinder=False):
    if   isParallel(face.Surface.Axis,FreeCAD.Vector(1,0,0),tol.tor_angle) :
        axis = FreeCAD.Vector(1,0,0)
    elif isParallel(face.Surface.Axis,FreeCAD.Vector(0,1,0),tol.tor_angle) :
        axis = FreeCAD.Vector(0,1,0)
    elif isParallel(face.Surface.Axis,FreeCAD.Vector(0,0,1),tol.tor_angle) :
        axis = FreeCAD.Vector(0,0,1)
        
    p1 = face.valueAt(0.,Vparams[0]) - face.Surface.Center
    z1 = p1.dot(axis)
    d1 = p1.cross(axis).Length
    
    p2 = face.valueAt(0.,Vparams[1])) - face.Surface.Center
    z2 = p2.dot(axis)
    d2 = p2.cross(axis).Length

    if isSameValue(z1,z2,tol.distance) :
       surfType = 'Plane'
       center   = face.Surface.Center + z1 * axis
       Vmid = (Vparams[0]+Vparams[1])*0.5
       pMid = face.valueAt(0,Vmid) - face.Surface.Center
       if pMid.dot(axis) < z1 :
           inSurf = True
       else:
           inSurf = False
       return (center,axis,face.Surface.MajorRadius,face.Surface.MajorRadius),surfType,inSurf        

    elif isSameValue(d1,d2,tol.distance) or forceCylinder :
       surfType = 'Cylinder'
       radius   = min(d1,d2)
       center   = face.Surface.Center
       if isSameValue(d1,face.Surface.MajorRadius,tol.distance) :
           Vmid = (Vparams[0]+Vparams[1])*0.5
           pMid = face.valueAt(0,Vmid) - center
           if pMid.cross(axis).Length < face.Surface.MajorRadius :
               inSurf = True
           else:
               inSurf = False
       else:
           if d1 < face.Surface.MajorRadius :
               inSurf = True
           else:    
               inSurf = False
       return (center,axis,radius,face.Surface.MinorRadius),surfType,inSurf
    
    else :
       surfType = 'Cone'
       za =  (z2*d1 - z1*d2)/(d1-d2)
       Apex   = face.Surface.Center + za * axis
       semiAngle = abs(math.atan(d1/(z1-za)))
       ConeAxis = axis if za < 0 else -axis 

       Vmid = (Vparams[0]+Vparams[1])*0.5
       pMid = face.valueAt(0,Vmid) - face.Surface.Center
       zMid = pMid.dot(axis) 
       dMid = pMid.cross(axis).Length
       
       dCone= d1 * (zMid-za)/(z1-za)
       inSurf = True if dMid < dCone  else False          
       return (Apex,ConeAxis,semiAngle,face.Surface.MinorRadius,face.Surface.MajorRadius),surfType,inSurf      



def cellDef(metaObj,Surfaces,UniverseBox):
  
    solids=metaObj.Solids
    delList = []

    PieceDef = BoolSequence(operator='OR')
    PieceObj = []
    for isol,solid in enumerate(solids):
        SurfPiece = []
        SurfObj   = []
        extraPlaneReverse = dict()

        flag_inv=isInverted(solid)
        solid_GU=GU.solid_GU(solid)
        lastTorus = -1
        for iface,face in enumerate(solid_GU.Faces):
            if abs(face.Area) < 1e-2 : continue 
            if face.Area < 0 :
               if opt.verbose : print('Warning : Negative surface Area')
            if face.Orientation not in ('Forward','Reversed') : continue
            if (flag_inv):
                orient_temp=face.Orientation
                if orient_temp=='Forward':
                    orient='Reversed'
                elif orient_temp=='Reversed':
                    orient='Forward'
            else:
                orient=face.Orientation

            surfaceType = str(face.Surface)
            if 'Sphere' in surfaceType : surfaceType = 'Sphere' 
            
            if (  surfaceType in ('<Cylinder object>','<Cone object>','Sphere') and orient=='Reversed'):
                if opt.verbose : print('Warning: Convex Cylinder, Sphere or Cone additional plane could be created')
                if (surfaceType == '<Cone object>'):
                                                                                      
                    if (isParallel(face.Surface.Axis,FreeCAD.Vector(1,0,0),tol.angle) or \
                        isParallel(face.Surface.Axis,FreeCAD.Vector(0,1,0),tol.angle) or \
                        isParallel(face.Surface.Axis,FreeCAD.Vector(0,0,1),tol.angle)):
                        if opt.verbose : print('Cone parallel to X, Y or Z not GQ is used')
                        id=getId(face.Surface,Surfaces) 
                        if (not(str(id) in SurfPiece)):
                            SurfPiece.append(str(id))
                            SurfObj.append(face)
                    else:
                        dim1 = face.Surface.Radius
                        dim2 = face.Surface.Radius
                        plane = GEOUNED_Surface(('Plane',(face.Surface.Apex,face.Surface.Axis,dim1,dim2)),UniverseBox,Face='Build')
                        id2,exist = Surfaces.addPlane(plane)
                        id = getId(face.Surface,Surfaces)

                        if exist :
                           p = Surfaces.getSurface(id2)
                           if isOposite(plane.Surf.Axis,p.Surf.Axis,tol.pln_angle):
                             id2= -id2
                        var='(%i:%i)' %(id,-id2)
                        if (not(var in SurfPiece)):
                            SurfPiece.append(var)
                            SurfObj.append([face,plane.shape])
                else:
                    id=getId(face.Surface,Surfaces)
                    if (not(str(id) in SurfPiece)):
                        SurfPiece.append(str(id))
                        SurfObj.append(face)
                        
                try:
                    plane=GenPlane(face,solid_GU)
                    if plane is not None:
                       plane=GU.Plane_GU(plane)    
                except:
                    plane=None
                    if opt.verbose : print('Warning: generation of additional plane has failed')
                
                idFace = id
                
                if plane is not None : 
                  p = GEOUNED_Surface(('Plane',(plane.Position,plane.Axis,plane.dim1,plane.dim2)),UniverseBox,Face='Build')
                  
                  id,exist = Surfaces.addPlane(p)
                  sign = signPlane(face.CenterOfMass,p)
                  if exist :
                     pp   = Surfaces.getSurface(id)
                     if isOposite(p.Surf.Axis,pp.Surf.Axis,tol.angle):
                          id= -id
                  id *= sign
                   
                  if ( idFace not in extraPlaneReverse.keys()):
                     extraPlaneReverse[idFace] = [str(id)]
                     SurfObj.append(p.shape)
                  else:
                     if str(id) not in extraPlaneReverse[idFace] :  
                       extraPlaneReverse[idFace].append(str(id))
                       SurfObj.append(p.shape)

                  #if (not(str(id) in SurfPiece)):
                     #SurfPiece.append(str(id))
                     #SurfObj.append(p.shape)

                if surfaceType == '<Cylinder object>' and False:    # temporary disabled doesn't work properly
                   PlanesId = ExtraPlaneCylFace(face,UniverseBox,Surfaces)
                   for id in PlanesId:
                     var = id
                     if (not(str(var) in SurfPiece)):
                       SurfPiece.append(str(var))
                       SurfObj.append(Surfaces.getSurface(abs(id)).shape)

            elif (surfaceType == '<Cone object>'): # this only if GQ is used
                if (isParallel(face.Surface.Axis,FreeCAD.Vector(1,0,0),tol.angle) or \
                    isParallel(face.Surface.Axis,FreeCAD.Vector(0,1,0),tol.angle) or \
                    isParallel(face.Surface.Axis,FreeCAD.Vector(0,0,1),tol.angle)):
                    #print 'Cone parallel to X, Y or Z not GQ is used'
                    id=getId(face.Surface,Surfaces)
                    id= -id 
                    if (not(str(id) in SurfPiece)):
                        SurfPiece.append(str(id))
                        SurfObj.append(face)
                else:
                    dim1 = face.Surface.Radius
                    dim2 = face.Surface.Radius
                    plane = GEOUNED_Surface(('Plane',(face.Surface.Apex,face.Surface.Axis,dim1,dim2)),UniverseBox,Face='Build')
                    id2,exist = Surfaces.addPlane(plane)
                    id1 = getId(face.Surface,Surfaces)

                    if exist :
                       p = Surfaces.getSurface(id2)
                       if isOposite(plane.Surf.Axis,p.Surf.Axis,tol.pln_angle):
                         id2= -id2
                    var='%i %i' %(-id1,id2)
                    if (not(var in SurfPiece)):
                        SurfPiece.append(var)
                        SurfObj.append([face,plane.shape])
                   
            elif (surfaceType == '<Toroid object>') :
                
                if (isParallel(face.Surface.Axis,FreeCAD.Vector(1,0,0),tol.angle) or \
                    isParallel(face.Surface.Axis,FreeCAD.Vector(0,1,0),tol.angle) or \
                    isParallel(face.Surface.Axis,FreeCAD.Vector(0,0,1),tol.angle)):
                 
                    id1=getId(face.Surface,Surfaces)
                    
                    if (orient == 'Forward'):

                        index,Uparams = solid_GU.TorusUParams[iface]
                        if index == lastTorus : continue
                        lastTorus = index
                        var = '-%i' %id1

                       # Old version seems not work porperly

                       # index,Uparams = solid_GU.TorusUParams[iface]
                       # if index == lastTorus : continue
                       # lastTorus = index
                       # closed,UminMax = Uparams
                       # if closed :
                       #    var = '-%i' %id1
                       # else:
                       #    (params1,params2),addOR = GenTorusAnnexUPlanes(face,UminMax)
                       #    plane = GEOUNED_Surface(('Plane',params1),UniverseBox,Face='Build')
                       #    ip1,exist = Surfaces.addPlane(plane)
                       #    #plane.shape.exportStep('p{}.stp'.format(ip1))
                       #    if exist :
                       #       p = Surfaces.getSurface(ip1)
                       #       if isOposite(plane.Surf.Axis,p.Surf.Axis,tol.pln_angle):
                       #          ip1=-ip1

                       #    if params2 : 
                       #       plane = GEOUNED_Surface(('Plane',params2),UniverseBox,Face='Build')
                       #       ip2,exist = Surfaces.addPlane(plane)
                       #       #plane.shape.exportStep('p{}.stp'.format(ip2))
                       #       if exist :
                       #          p = Surfaces.getSurface(ip2)
                       #          if isOposite(plane.Surf.Axis,p.Surf.Axis,tol.pln_angle):
                       #             ip2=-ip2

                       #       if addOR :                      
                       #          var = '-%i (%i:%i)' %(id1,ip1,ip2)
                       #       else:   
                       #          var = '-%i %i %i' %(id1,ip1,ip2)
                       #    else:
                       #        var = '-%i %i' %(id1,ip1)
                    else:
                        index,Vparams = solid_GU.TorusVParams[iface]
                        if index == lastTorus : continue
                        lastTorus = index
                        closed,VminMax = Vparams
                        if closed :
                           var = '%i' %id1
                        else:
                           surfParams,surfType,inSurf = GenTorusAnnexVSurface(face,VminMax,opt.forceCylinder)

                           if surfType == 'Cone' :
                              cone = GEOUNED_Surface(('Cone',surfParams),UniverseBox,Face='Build')
                              id2,exist = Surfaces.addCone(cone)
                                                     
                           elif surfType == 'Cylinder' :
                              cyl = GEOUNED_Surface(('Cylinder',surfParams),UniverseBox,Face='Build')
                              id2,exist = Surfaces.addCylinder(cyl)
                              #cyl.shape.exportStep('cyl{}.stp'.format(id2))
                                                    
                           elif surfType == 'Plane' :
                              plane = GEOUNED_Surface(('Plane',surfParams),UniverseBox,Face='Build')
                              id2,exist = Surfaces.addPlane(plane)
                              if exist :
                                 p = Surfaces.getSurface(id2)
                                 if isOposite(plane.Surf.Axis,p.Surf.Axis,tol.pln_angle):
                                     id2=-id2

                           var = '%i %i' %(id1,-id2 if inSurf else id2)

                    if (not(var in SurfPiece)):
                        SurfPiece.append(var)
                        SurfObj.append(face)
                else:
                    if opt.verbose : print('Only Torus with axis along X, Y , Z axis can be reproduced')
            else:
                id=getId(face.Surface,Surfaces)
                surf = face
                if id == 0:
                    if opt.verbose : print('Warning: ', surfaceType, ' not found in surface list')
                    if surfaceType == '<Plane object>' :
                      dim1 = face.ParameterRange[1]-face.ParameterRange[0]
                      dim2 = face.ParameterRange[3]-face.ParameterRange[2]  
                      plane = GEOUNED_Surface(('Plane',(face.Surface.Position,face.Surface.Axis,dim1,dim2)),UniverseBox,Face='Build')
                      id,exist    = Surfaces.addPlane(plane)
                      surf  = plane.shape
                    elif surfaceType == '<Cylinder object>':
                      dimL = face.ParameterRange[3]-face.ParameterRange[2]  
                      cylinder = GEOUNED_Surface(('Cylinder',(face.Surface.Center,face.Surface.Axis,face.Surface.Radius,dimL)),UniverseBox,Face='Build')
                      id,exist    = Surfaces.addCylinder(cylinder)
                      surf  = cylinder.shape
               
                if (orient == 'Reversed'):
                    var=id
                elif (orient == 'Forward'):
                    var=-id

                if  surfaceType == '<Plane object>':
                   s = Surfaces.getSurface(id)
                   if  isOposite(face.Surface.Axis,s.Surf.Axis,tol.pln_angle) :
                      var=-var
                       
                if (str(var) in SurfPiece):
                    continue

                if (not(str(var) in SurfPiece)):
                    SurfPiece.append(str(var))
                    SurfObj.append(surf)

        if extraPlaneReverse :
          for extra in extraPlaneReverse.values() :
             if len(extra) ==1 :
                if extra[0] not in SurfPiece:
                   SurfPiece.append(extra[0])
             else:
                SurfPiece.append('({})'.format(':'.join(extra)))

        SurfPieceBool = BoolSequence(' '.join(SurfPiece))
        # possible expresion for e
        #  i1
        #  i1 i2SurfPiece
        #  i1 i2 i3
        #  (i1:i2)
        #  i1 (i2:i3)
        
 #         semi = e.find(':')
 #         blk  = e.find(' ')
 #         print (e)
 #         if semi != -1 :
 #           orTerm = expOR(int(e[1:semi]),int(e[semi+1:-1]))
 #           SurfPieceBool.add(orTerm)  
 #         elif blk != -1 :  
 #           SurfPieceBool.add(int(e[1:blk]),int(e[blk+1:-1]))
 #         else :
 #          SurfPieceBool.add(int(e))

        if SurfPieceBool.elements : 
          PieceDef.append(SurfPieceBool)
          PieceObj.append(SurfObj)
        else:
          delList.append(isol)

    for isol in reversed(delList) :
       del metaObj.Solids[isol]
    metaObj.setDefinition(PieceDef)
    metaObj.setFaces(PieceObj)
    return 
    
def ExtraPlaneCylFace(face,Box,Surfaces):
    wire = face.OuterWire
    planesId = []
    for e in wire.OrderedEdges :
       curve = str(e.Curve)
       if (curve[0:6] == 'Circle' or curve == '<Ellipse object>' ):
          dir = e.Curve.Axis
          center = e.Curve.Center
          if curve == '<Ellipse object>' :
            dim1   = e.Curve.MinorRadius
            dim2   = e.Curve.MajorRadius
          else:    
            dim1   = e.Curve.Radius
            dim2   = e.Curve.Radius
          plane = GEOUNED_Surface(('Plane',(center,dir,dim1,dim2)),Box,Face='Build')
          sign = signPlane(face.CenterOfMass,plane)
          id,exist = Surfaces.addPlane(plane)
          if exist :
             pp   = Surfaces.getSurface(id)
             if isOposite(plane.Surf.Axis,pp.Surf.Axis,tol.pln_angle):
                 id=-id
          planesId.append(id)
    return planesId
         
