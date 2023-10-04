#
# Set of useful functions used in different parts of the code
#
import math

def isSameValue(v1,v2,tolerance=1e-6):
    return abs(v1-v2) < tolerance 


def isOposite(Vector1,Vector2,tolerance=1e-6):
    return abs(Vector1.getAngle(-Vector2)) < tolerance
        

def isParallel(Vector1,Vector2,tolerance=1e-6):
    angle = abs(Vector1.getAngle(Vector2))
    return (angle < tolerance or isSameValue(angle,math.pi,tolerance) ) 


def isInLine(Point,Dir,Pnt_line,tolerance=1e-6):
    r12 = Point - Pnt_line
    return (isParallel(Dir,r12) or (r12.Length < tolerance))      

    
def isInPoints(point,Points,tolerance=1e-5):
   if len(Points)>0:
        for P in Points:
            if point.isEqual(P,tolerance):
                return True
   return False

def isInEdge(Edge1,Edge2,tolerance=1e-8):
    Ver1=Edge1.Vertexes
    Ver2=Edge2.Vertexes
    con1 = Ver1[0].Point.isEqual(Ver2[0].Point,tolerance) or Ver1[0].Point.isEqual(Ver2[1].Point,tolerance)
    con2 = Ver1[1].Point.isEqual(Ver2[0].Point,tolerance) or Ver1[1].Point.isEqual(Ver2[1].Point,tolerance)    
    return  (con1 and con2)


def isInPlane(Point,plane,dtolerance=1e-7):
    return abs(Point.distanceToPlane(plane.Surf.Position,plane.Surf.Axis)) < dtolerance


def isInTolerance(val,tol,fuzzyLow,fuzzyHigh):
    if abs(val) < fuzzyLow :
        return True,False   # 1) isintolerance 2) fuzzy
    elif abs(val) < tol :
        return True,True
    elif abs(val) > fuzzyHigh :
        return False,False
    else:
        return False,True

def signPlane(Point,plane):
    value = plane.Surf.Axis.dot(Point-plane.Surf.Position) 
    if (value>=0.0):
        sign=1
    else:
        sign=-1
    return sign

class PlaneParams:
   def __init__(self,params,real=True):
      self.Position = params[0] 
      self.Axis     = params[1]
      self.dimL     = params[2]
      self.solidID  = params[3]
      self.real     = real
      

   def __str__(self):
      outstr = '''Plane :
    Axis     : {}  {}  {} 
    Position : {}  {}  {}'''.format(self.Axis.x    ,self.Axis.y    ,self.Axis.z,\
                                    self.Position.x,self.Position.y,self.Position.z)
      return outstr

class CylinderParams:
   def __init__(self,params,real=True):
      self.Center = params[0] 
      self.Axis   = params[1]
      self.Radius = params[2]
      self.dimL   = params[3]
      self.real   = real

   def __str__(self):
      outstr = '''Cylinder :
    Axis     : {}  {}  {} 
    Center   : {}  {}  {}
    Radius   : {}  '''.format(self.Axis.x    ,self.Axis.y    ,self.Axis.z,\
                             self.Center.x,self.Center.y,self.Center.z,\
                              self.Radius)
      return outstr

class ConeParams:
   def __init__(self,params,real=True):
      self.Apex   = params[0] 
      self.Axis   = params[1]
      self.SemiAngle = params[2]
      self.dimL   = params[3]
      self.dimR   = params[4]
      self.real   = real

   def __str__(self):
      outstr = '''Cone :
    Axis     : {}  {}  {} 
    Center   : {}  {}  {}
    SemiAngle: {}  '''.format(self.Axis.x    ,self.Axis.y    ,self.Axis.z,\
                             self.Apex.x,self.Apex.y,self.Apex.z,\
                              self.SemiAngle)
      return outstr


class SphereParams:
   def __init__(self,params):
      self.Center = params[0] 
      self.Radius = params[1]

   def __str__(self):
      outstr = '''Sphere :
    Center   : {}  {}  {}
    Radius   : {}  '''.format(self.Center.x,self.Center.y,self.Center.z,\
                              self.Radius)
      return outstr

class TorusParams:
    def __init__(self,params):
      self.Center      = params[0]
      self.Axis        = params[1]
      self.MajorRadius = params[2]
      self.MinorRadius = params[3]

    def __str__(self):
      outstr = '''Torus :
    Axis     : {}  {}  {} 
    Center   : {}  {}  {}
    MajorRadius: {}
    MinorRadius: {} '''.format(self.Axis.x    ,self.Axis.y    ,self.Axis.z,\
                              self.Center.x,self.Center.y,self.Center.z,\
                              self.MajorRadius,self.MinorRadius)
      return outstr
       
