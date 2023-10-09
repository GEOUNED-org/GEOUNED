#Init file of GEOUNED module
#

# We load the STEP and the materials
import sys
import FreeCAD,Part
import configparser

from os import path,mkdir
from datetime import datetime


import GEOUNED.LoadFile.LoadSTEP as Load
import GEOUNED.Decompose.Decom_one as Decom
import GEOUNED.Utils.Functions as UF
import GEOUNED.Conversion.CellDefinition as Conv
from GEOUNED.Write.Functions import writeMCNPCellDef
from GEOUNED.Write.WriteFiles import writeGeometry
from GEOUNED.CodeVersion import *
from GEOUNED.Utils.Options.Classes import Tolerances, MCNP_numeric_format, Options 
from GEOUNED.Utils.BooleanSolids import buildCTableFromSolids
from GEOUNED.Cuboid.translate import translate

import GEOUNED.Void.Void as Void


class GEOUNED() :

   def __init__(self,title='Geouned conversion'):
      """ """
      self.__dict__ = dict() 
      self.__dict__['stepFile']     = ''
      self.__dict__['geometryName'] = ''
      self.__dict__['matFile']      = ''
      self.__dict__['outFormat']    = ('mcnp',)
      self.__dict__['title']        = title
      self.__dict__['voidGen']      = True
      self.__dict__['debug']        = False
      self.__dict__['compSolids']   = True
      self.__dict__['volSDEF']      = False
      self.__dict__['dummyMat']     = False
      self.__dict__['volCARD']      = True
      self.__dict__['UCARD']        = None
      self.__dict__['simplify']     = 'No'
      self.__dict__['cellRange']    = []
      self.__dict__['exportSolids'] = ''
      self.__dict__['minVoidSize']  = 200  # units mm
      self.__dict__['maxSurf']      = 50
      self.__dict__['maxBracket']   = 30
      self.__dict__['voidMat']      = []
      self.__dict__['voidExclude']  = []
      self.__dict__['startCell']    = 1
      self.__dict__['startSurf']    = 1
      self.__dict__['cellCommentFile'] = False
      self.__dict__['cellSummaryFile'] = True
      self.__dict__['sortEnclosure'] = False

   def SetOptions(self):
      toleranceKwrd = ( 'relativeTolerance', 'relativePrecision', 'singleValue', 'generalDistance', 'generalAngle',
                        'planeDistance', 'planeAngle', 'cylinderDistance', 'cylinderAngle', 'sphereDistance', 'coneDistance',
                        'coneAngle', 'torusDistance', 'torusAngle')
      numericKwrd = ('P_abc', 'P_d', 'P_xyz', 'S_r', 'S_xyz', 'C_r', 'C_xyz', 'K_tan2', 'K_xyz', 'T_r', 'T_xyz', 'GQ_1to6', 'GQ_7to9', 'GQ_10' )
      tolKwrdEquiv = {'relativeTolerance':'relativeTol' , 'relativePrecision':'relativePrecision', 'singleValue':'value', 'generalDistance':'distance',
                      'generalAngle': 'angle', 'planeDistance':'pln_distance', 'planeAngle':'pln_angle', 'cylinderDistance':'cyl_distance', 
                      'cylinderAngle':'cyl_angle', 'sphereDistance':'sph_distance', 'coneDistance':'kne_distance', 'coneAngle': 'kne_angle',
                      'torusDistance':'tor_distance', 'torusAngle':'tor_angle'}

      config = configparser.ConfigParser()
      config.optionxform = str
      config.read(self.__dict__['title'])
      for section in config.sections():
         if section == 'Files':
            for key in config['Files'].keys() :
                 if key in ('geometryName','matFile','title') :
                     self.set(key, config.get('Files',key))   

                 elif key == 'stepFile':
                     value = config.get('Files',key).strip()
                     if value[0] in ('(','[')  and value[-1] in (']',')') :
                        data=value[1:-1].split(',')
                        self.set(key,data)
                     else :
                        self.set(key, value)   

                 elif key == 'outFormat':
                     raw = config.get('Files',key).strip()
                     values = tuple( x.strip() for x in raw.split(',') )
                     outFormat = []
                     for v in values:
                        if   v.lower() == 'mcnp'       : outFormat.append('mcnp')
                        elif v.lower() == 'openmc_xml' : outFormat.append('openMC_XML')
                        elif v.lower() == 'openmc_py'  : outFormat.append('openMC_PY')
                     self.set(key, tuple(outFormat))   
                     

         elif section == 'Parameters':
            for key in config['Parameters'].keys() :
                 if key in ('voidGen','debug','compSolids','volSDEF','volCARD',
                            'dummyMat','cellSummaryFile','cellCommentFile','sortEnclosure') :
                     self.set(key, config.getboolean('Parameters',key))   
                 elif key in ('minVoidSize','maxSurf','maxBracket','startCell','startSurf') :
                     self.set(key, config.getint('Parameters',key))   
                 elif key in ('exportSolids','UCARD','simplify') :
                     self.set(key, config.get('Parameters',key))   
                 elif key == 'voidMat':
                     value = config.get('Parameters',key).strip()
                     data=value[1:-1].split(',')
                     self.set(key,(int(data[0]),float(data[1]),data[2]))
                 else:
                     value = config.get('Parameters',key).strip()
                     data=value[1:-1].split(',')
                     self.set(key,tuple(map(int,data)))

         elif section == 'Options':
            for key in config['Options'].keys() :
               if key in ('forceCylinder','newSplitPlane','delLastNumber','verbose','quadricPY','Facets','prnt3PPlane') :
                   setattr(Options,key, config.getboolean('Options',key))
               elif key in ('enlargeBox','nPlaneReverse','splitTolerance'):
                   setattr(Options,key, config.getfloat('Options',key))
         
         elif section == 'Tolerances':
            for key in config['Tolerances'].keys() :
                 if key == 'relativeTolerance' :
                     setattr(Tolerances,key, config.getboolean('Tolerances',key))
                 elif key in toleranceKwrd : setattr(Tolerances,tolKwrdEquiv[key], config.getfloat('Tolerances',key))   

         elif section == 'MCNP_Numeric_Format':
            for key in config['MCNP_Numeric_Format'].keys() :
                 if key in numericKwrd : setattr(MCNP_numeric_format,key, config.get('MCNP_Numeric_Format',key))   
             
         else:
            print('bad section name : {}'.format(section))

      if self.__dict__['geometryName'] == '' :
         self.__dict__['geometryName'] = self.__dict__['stepFile'][:-4]
      print(self.__dict__)

   def set(self,kwrd,value):

       if kwrd not in self.__dict__.keys():
          print('Bad entry : {}'.format(kwrd))
          return

       if kwrd == 'stepFile' :
          if isinstance(value,(list,tuple)) :
             for v in value :
               if not isinstance(v,str) :
                  print('elemt in {} list should be string'.format(kwrd))
                  return
          elif not isinstance(value,str) :
                print('{} should be string or tuple of strings'.format(kwrd))
                return

       elif kwrd == 'UCARD' :
           if value == 'None' :
              value = None
           elif value.isdigit() :
              value = int(value)
           else:
              print('{} value should be None or integer'.format(kwrd))
              return
       elif kwrd == 'outFormat' :
           if len(value) == 0: return
       elif kwrd in ('geometryName','matFile','exportSolids'):
           if not isinstance(value,str) :
              print('{} value should be str instance'.format(kwrd))
              return   
       elif kwrd in ('cellRange','voidMat','voidExclude'):
           if not isinstance(value,(list,tuple)) :
              print('{} value should be list or tuple'.format(kwrd))
              return
       elif kwrd in ('minVoidSize','maxSurf','maxBracket','startCell','startSurf'):
           if not isinstance(value,int) :
              print('{} value should be integer'.format(kwrd))
              return
       elif kwrd in ('voidGen','debug','compSolids','simplifyCTable','volSDEF','volCARD', 
                     'dummyMat','cellSummaryFile','cellCommentFile','sortEnclosure') :
           if not isinstance(value,bool) :
              print('{} value should be boolean'.format(kwrd))
              return    

       self.__dict__[kwrd] = value
       if kwrd == 'stepFile' and  self.__dict__['geometryName'] == '' :
            if isinstance(value,(tuple,list)):
               self.__dict__['geometryName'] == 'joined_step_files' 
            else: 
               self.__dict__['geometryName'] == value[:-4]

             
   def Start(self):

       print('start')
       FreeCAD_Version = '{V[0]:}.{V[1]:}.{V[2]:}'.format(V=FreeCAD.Version())
       print('GEOUNED version {} {} \nFreeCAD version {}'.format(GEOUNED_Version,GEOUNED_ReleaseDate,FreeCAD_Version))

       code_setting = self.__dict__
       if code_setting is None:
          print('Cannot run the code. Input are missing')
          exit()
       if code_setting['stepFile'] == '':
          print('Cannot run the code. Step file name is missing')
          exit()
          
       stepfile = code_setting['stepFile']
       matfile  = code_setting['matFile']
       startTime = datetime.now()

       if isinstance(stepfile,(list,tuple)):
          MetaChunk      = []
          EnclosureChunk = []
          for stp in stepfile :
             print('read step file : {}'.format(stp))   
             Meta,Enclosure = Load.LoadCAD(stp,matfile)
             MetaChunk.append(Meta)
             EnclosureChunk.append(Enclosure)
          MetaList = joinMetaLists(MetaChunk)
          EnclosureList = joinMetaLists(EnclosureChunk)
       else:
          print('read step file : {}'.format(stepfile)) 
          MetaList,EnclosureList = Load.LoadCAD(stepfile,matfile,code_setting['voidMat'],code_setting['compSolids'])

       print('End of loading phase')
       tempstr1 = str(datetime.now()-startTime)
       print(tempstr1)
       tempTime = datetime.now()

       # Select a specific solid range from original STEP solids
       if code_setting['cellRange'] :
           MetaList = MetaList[code_setting['cellRange'][0]:code_setting['cellRange'][1]]

       # export in STEP format solids read from input file
       # terminate excution
       if code_setting['exportSolids'] != '':    
          solids=[]
          for m in MetaList:
              if m.IsEnclosure: continue 
              solids.extend(m.Solids)
          Part.makeCompound(solids).exportStep(code_setting['exportSolids'])
          print ('Solids exported in file :{}'.format(code_setting['exportSolids']))
          print ('GEOUNED Finish. No solid translation performed.')
          sys.exit()

       # set up Universe
       if EnclosureList :
          UniverseBox = getUniverse(MetaList+EnclosureList)
       else:
          UniverseBox = getUniverse(MetaList)
       Comsolids= []

       surfOffset = code_setting['startSurf'] - 1
       Surfaces    = UF.Surfaces_dict(offset=surfOffset)


       warnSolids     = []
       warnEnclosures = []
       tempTime0 = datetime.now()
       if not Options.Facets:

          # decompose all solids in elementary solids (convex ones)
          warningSolidList = DecomposeSolids(MetaList,Surfaces,UniverseBox,code_setting,True)

          # decompose Enclosure solids
          if code_setting['voidGen'] and EnclosureList :
             warningEnclosureList = DecomposeSolids(EnclosureList,Surfaces,UniverseBox,code_setting,False)
           
          print('End of decomposition phase')

          # start Building CGS cells phase   
 
          for j,m in enumerate(MetaList):
              if m.IsEnclosure : continue
              print('Building cell: ', j)
              Conv.cellDef(m,Surfaces,UniverseBox)
              if j in warningSolidList:
                  warnSolids.append(m)
              if not m.Solids:
                 print('none',j,m.__id__)
                 print(m.Definition)

       else:
          translate(MetaList,Surfaces,UniverseBox,code_setting)
          # decompose Enclosure solids
          if code_setting['voidGen'] and EnclosureList :
             warningEnclosureList = DecomposeSolids(EnclosureList,Surfaces,UniverseBox,code_setting,False)

       tempstr2 = str(datetime.now()-tempTime)
       print(tempstr2)
         

       #  building enclosure solids       

       if code_setting['voidGen'] and EnclosureList :
           for j,m in enumerate(EnclosureList):
             print('Building Enclosure Cell: ', j)
             Conv.cellDef(m,Surfaces,UniverseBox)
             if j in warningEnclosureList:
                warnEnclosures.append(m)
        

       tempTime1 = datetime.now()

       printWarningSolids(warnSolids,warnEnclosures)

       # void generation phase
       MetaVoid = []
       if code_setting['voidGen'] :
         print('Build Void')
         print(code_setting['voidExclude'])
         if not code_setting['voidExclude']:
            MetaReduced = MetaList            
         else :
            MetaReduced = excludeCells(MetaList,code_setting['voidExclude'])

         if MetaList :
            init = MetaList[-1].__id__ -len(EnclosureList)
         else:
            init = 0
         MetaVoid = Void.voidGeneration(MetaReduced,EnclosureList,Surfaces,UniverseBox,code_setting,init)

       if code_setting['simplify'] == 'full' :
           Surfs = {}
           for lst in Surfaces.values():
              for s in lst:
                Surfs[s.Index] = s

           for c in MetaList:
              if c.Definition.level == 0 or c.IsEnclosure: continue
              print('simplify cell',c.__id__)
              Box =  getBox(c)
              CT = buildCTableFromSolids(Box,(c.Surfaces,Surfs),option='full')
 
                   
              c.Definition.simplify(CT)
              res = c.Definition.clean()
              if res is not None:
                  print('unexpected constant cell {} :{}'.format(cell.__id__,res)) 



       tempTime2 = datetime.now()
       print('build Time:',tempTime2-tempTime1)
 
       print(datetime.now()-startTime)

       cellOffSet = code_setting['startCell'] - 1
       if EnclosureList and code_setting['sortEnclosure'] :
          # sort group solid cell / void cell sequence in each for each enclosure
          # if a solid belong to several enclosure, its definition will be written 
          # for the highest enclosure level or if same enclosure level in the first 
          # enclosure found
          MetaList = sortEnclosure(MetaList,MetaVoid,cellOffSet) 
       else:
          # remove Null Cell and apply cell numbering offset
          deleted = []
          icount = 0
          for i,m in enumerate(MetaList):
              if m.NullCell or m.IsEnclosure :
                deleted.append(i)
                icount += 1
                continue
              m.__id__ += cellOffSet - icount
          for i in reversed(deleted):
             del MetaList[i]       

          lineComment = """\
##########################################################
             VOIDS 
##########################################################"""
          mc = UF.GEOUNED_Solid(None)
          mc.Comments = lineComment
          MetaList.append(mc)

          deleted = []
          icount = 0
          for i,m in enumerate(MetaVoid):
              if m.NullCell :
                deleted.append(i)
                icount += 1
                continue
              m.__id__ += cellOffSet - icount
          for i in reversed(deleted):
             del MetaVoid[i]       

          MetaList.extend(MetaVoid)

       # write outputformat input
       writeGeometry(UniverseBox,MetaList,Surfaces,code_setting)
   
       print('End of MCNP translation phase')

       print('Process finished')
       print(datetime.now()-startTime)
    
       print('Translation time of solid cells', tempTime1-tempTime0)
       print('Translation time of void cells', tempTime2-tempTime1)


def DecomposeSolids(MetaList,Surfaces,UniverseBox,setting,meta):
    totsolid = len(MetaList)
    warningSolids = []
    for i,m in enumerate(MetaList):
        if meta and m.IsEnclosure: continue
        print('Decomposing solid: {}/{} '.format(i,totsolid))
        if setting['debug'] :
            print(m.Comments)
            if m.IsEnclosure :
               m.Solids[0].exportStep('origEnclosure_{}.stp'.format(i))  
            else:    
               m.Solids[0].exportStep('origSolid_{}.stp'.format(i))
        
        comsolid,err=Decom.SplitSolid(Part.makeCompound(m.Solids),UniverseBox)

        if err != 0 :
           if not path.exists('Suspicious_solids') :
              mkdir("Suspicious_solids")
           if m.IsEnclosure:   
              m.Solids[0].exportStep(u'Suspicious_solids/Enclosure_original_{}.stp'.format(i))
              comsolid.exportStep(   u'Suspicious_solids/Enclosure_split_{}.stp'.format(i))
           else: 
              m.Solids[0].exportStep(u'Suspicious_solids/Solid_original_{}.stp'.format(i))
              comsolid.exportStep(   u'Suspicious_solids/Solid_split_{}.stp'.format(i))

           warningSolids.append(i)
           
        if setting['debug'] :
            if m.IsEnclosure :
                comsolid.exportStep('compEnclosure_{}.stp'.format(i))
            else:
                comsolid.exportStep('compSolid_{}.stp'.format(i))
        Surfaces.extend(Decom.ExtractSurfaces(comsolid,'All',UniverseBox,MakeObj=True))
        m.setCADSolid()
        m.updateSolids(comsolid.Solids)

    return warningSolids   


def getBox(comp):
    bb=FreeCAD.BoundBox(comp.BoundBox)
    apexList = findCone(comp.Faces)
    if apexList :
       xMin,xMax,yMin,yMax,zMin,zMax = bb.XMin, bb.XMax, bb.YMin, bb.YMax, bb.ZMin, bb.ZMax
       for apex in apexList:
          if apex.x < xMin   : xMin = apex.x
          elif apex.x > xMax : xMax = apex.x
          if apex.y < yMin   : yMin = apex.y
          elif apex.y > yMax : yMax = apex.y
          if apex.z < zMin   : zMin = apex.z
          elif apex.z > zMax : zMax = apex.z
       xLength = xMax-xMin + Options.enlargeBox
       yLength = yMax-yMin + Options.enlargeBox
       zLength = zMax-zMin + Options.enlargeBox
       xMin -= 0.5*Options.enlargeBox
       yMin -= 0.5*Options.enlargeBox
       zMin -= 0.5*Options.enlargeBox
    else:
       bb.enlarge(Options.enlargeBox)
       xMin,yMin,zMin = bb.XMin, bb.YMin, bb.ZMin
       xLength,yLength,zLength = bb.XLength,  bb.YLength,  bb.ZLength 
       
    return Part.makeBox( xLength, yLength, zLength,
                         FreeCAD.Vector( xMin, yMin, zMin),
                         FreeCAD.Vector(0,0,1))

def findCone(Faces):
   coneObj = []
   for list1 in Faces:
     for f in list1:
        if type(f) is list :
           for ff in f :
              if ff is not None :
                 if "Cone" in str(ff.Surface) : coneObj.append(ff)
        else:
           if "Cone" in str(f.Surface) : coneObj.append(f)

   apexList = []
   for c in coneObj :
      apexList.append(c.Surface.Apex)
   return apexList

def getUniverse(MetaList):
      d = 10
      Box  = MetaList[0].BoundBox
      xmin = Box.XMin
      xmax = Box.XMax
      ymin = Box.YMin
      ymax = Box.YMax
      zmin = Box.ZMin
      zmax = Box.ZMax
      for m in  MetaList[1:]:
       # MIO. This was removed since in HELIAS the enclosure cell is the biggest one
       # if m.IsEnclosure: continue
        xmin = min(m.BoundBox.XMin,xmin)
        xmax = max(m.BoundBox.XMax,xmax)
        ymin = min(m.BoundBox.YMin,ymin)
        ymax = max(m.BoundBox.YMax,ymax)
        zmin = min(m.BoundBox.ZMin,zmin)
        zmax = max(m.BoundBox.ZMax,zmax)
     
      return FreeCAD.BoundBox(FreeCAD.Vector(xmin-d,ymin-d,zmin-d),FreeCAD.Vector(xmax+d,ymax+d,zmax+d))


def printWarningSolids(warnSolids,warnEnclosures):

    if warnSolids or warnEnclosures :
       fic = open('Warning_Solids_definition.txt','w')
    else:
       return

    if warnSolids :
       lines = 'Solids :\n'
       for sol in warnSolids:
          lines += '\n'
          lines += '{}\n'.format(sol.__id__)
          lines += '{}\n'.format(sol.Comments)
          lines += '{}\n'.format(writeMCNPCellDef(sol.Definition))
       fic.write(lines)

    if warnEnclosures :
       lines = 'Enclosures :\n'
       for sol in warnEnclosures:
          lines += '\n'
          lines += '{}\n'.format(sol.__id__)
          lines += '{}\n'.format(sol.Comments)
          lines += '{}\n'.format(writeMCNPCellDef(sol.Definition))

       fic.write(lines)

    fic.close()


def joinMetaLists(MList):

   newMetaList = MList[0]
   if MList[0]:
      for M in MList[1:] :
         lastID = newMetaList[-1].__id__+1
         for i,meta in enumerate(M) :
            meta.__id__ = lastID + i
            newMetaList.append(meta)

   return newMetaList     

def excludeCells(MetaList,labelList):
    voidMeta = []
    for m in MetaList:
        if m.IsEnclosure : continue
        found = False
        for label in labelList:
           if label in m.Comments :
              found = True
              break
        if not found :
           voidMeta.append(m)

    return voidMeta       



def sortEnclosure(MetaList,MetaVoid,offSet=0):

    newList = {}
    for m in MetaVoid:
       if m.EnclosureID in newList.keys():
          newList[m.EnclosureID].append(m)
       else:
          newList[m.EnclosureID] = [m]
    
    icount = offSet
    newMeta = []
    for m in MetaList:
      if m.NullCell : continue
      if m.IsEnclosure:
         lineComment = """\
##########################################################
             ENCLOSURE {}
##########################################################""".format(m.EnclosureID)
         mc = UF.GEOUNED_Solid(None)
         mc.Comments = lineComment
         newMeta.append(mc)
         for e in newList[m.EnclosureID] :
            if e.NullCell : continue
            icount += 1
            e.__id__ = icount
            newMeta.append(e)
         lineComment = """\
##########################################################
            END  ENCLOSURE {}
##########################################################""".format(m.EnclosureID)
         mc = UF.GEOUNED_Solid(None)
         mc.Comments = lineComment
         newMeta.append(mc)

      else:
         icount += 1
         m.__id__ = icount
         newMeta.append(m)


    lineComment = """\
##########################################################
             VOIDS 
##########################################################"""
    mc = UF.GEOUNED_Solid(None)
    mc.Comments = lineComment
    newMeta.append(mc)

    for v in newList[0] :
        if v.NullCell : continue
        icount += 1
        v.__id__ = icount
        newMeta.append(v)
    
    return newMeta

def sortEnclosure_old(MetaList):
    newList = {}
    enclID = []
    for m in MetaList:
      if m.ParentEnclosureID not in enclID:
         enclID.append(m.ParentEnclosureID)
         newList[m.ParentEnclosureID]=[m]
      else:
         newList[m.ParentEnclosureID].append(m)
    
    if -2 in enclID:
       enclID.remove(-2)
       graveYard = True
    else:
       graveYard = False

    icount = 0
    newMeta = []
    for eID in enclID:
      if eID != -1 :
         lineComment = """\
##########################################################
             ENCLOSURE {}
##########################################################""".format(eID)
         mc = UF.GEOUNED_Solid(None)
         mc.Comments = lineComment
         newMeta.append(mc)

      for m in newList[eID] :
         icount += 1
         m.__id__ = icount
         newMeta.append(m)

    if graveYard:
      for m in newList[-2] :
         icount += 1
         m.__id__ = icount
         newMeta.append(m)
    
    return newMeta





