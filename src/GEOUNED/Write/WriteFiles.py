import GEOUNED.Write.AdditionalFiles as OutFiles
from GEOUNED.Write.MCNPFormat        import MCNP_input
from GEOUNED.Write.OpenMCFormat      import OpenMC_input 

def writeGeometry(UniverseBox,MetaList,Surfaces,code_setting):

   baseName  = code_setting['geometryName']

   # write cells comments in file
   if code_setting['cellCommentFile']:
       OutFiles.commentsWrite(baseName,MetaList)
   if code_setting['cellSummaryFile']:
       OutFiles.summaryWrite(baseName,MetaList)

   
   if 'mcnp' in code_setting['outFormat']:
       mcnpFilename = baseName + '.mcnp'
       outBox = (UniverseBox.XMin,UniverseBox.XMax,UniverseBox.YMin,UniverseBox.YMax,UniverseBox.ZMin,UniverseBox.ZMax)
       if code_setting['voidGen'] :
          outSphere = (Surfaces['Sph'][-1].Index,Surfaces['Sph'][-1].Surf.Radius)
       else:
          outSphere = None       

       MCNPfile = MCNP_input(MetaList,Surfaces,code_setting)
       MCNPfile.setSDEF((outSphere,outBox))
       MCNPfile.writeInput(mcnpFilename)

   if 'openMC_XML' in code_setting['outFormat'] or \
      'openMC_PY'  in code_setting['outFormat'] :
       OMCFile = OpenMC_input(MetaList,Surfaces,code_setting)

   if 'openMC_XML' in code_setting['outFormat']:
       omcFilename = baseName + '.xml'
       OMCFile.writeXML(omcFilename)
       

   if 'openMC_PY' in code_setting['outFormat']:
       omcFilename = baseName + '.py'
       OMCFile.writePY(omcFilename)

   if 'serpent' in code_setting['outFormat']:
       mcnpFilename = baseName + '.serp'
       outBox = (UniverseBox.XMin,UniverseBox.XMax,UniverseBox.YMin,UniverseBox.YMax,UniverseBox.ZMin,UniverseBox.ZMax)
       if code_setting['voidGen'] :
          outSphere = (Surfaces['Sph'][-1].Index,Surfaces['Sph'][-1].Surf.Radius)
       else:
          outSphere = None       

       MCNPfile = MCNP_input(MetaList,Surfaces,code_setting)
       MCNPfile.setSDEF((outSphere,outBox))
       MCNPfile.writeInput(mcnpFilename)   
