from . import AdditionalFiles as OutFiles
from .MCNPFormat import MCNP_input
from .OpenMCFormat import OpenMC_input
from .PHITSFormat import PHITS_input
from .SerpentFormat import Serpent_input


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
       serpentFilename = baseName + '.serp'
       outBox = (UniverseBox.XMin,UniverseBox.XMax,UniverseBox.YMin,UniverseBox.YMax,UniverseBox.ZMin,UniverseBox.ZMax)
       if code_setting['voidGen'] :
          outSphere = (Surfaces['Sph'][-1].Index,Surfaces['Sph'][-1].Surf.Radius)
       else:
          outSphere = None       

       Serpentfile = Serpent_input(MetaList,Surfaces,code_setting)
       # Serpentfile.setSDEF((outSphere,outBox))
       Serpentfile.writeInput(serpentFilename)   

   if 'phits' in code_setting['outFormat']:
       phitsFilename = baseName + '.inp'
       PHITS_outBox = (UniverseBox.XMin,UniverseBox.XMax,UniverseBox.YMin,UniverseBox.YMax,UniverseBox.ZMin,UniverseBox.ZMax)
       if code_setting['voidGen'] :
          PHITS_outSphere = (Surfaces['Sph'][-1].Index,Surfaces['Sph'][-1].Surf.Radius)
       else:
          PHITS_outSphere = None       

       PHITSfile = PHITS_input(MetaList,Surfaces,code_setting)
       # PHITSfile.setSDEF_PHITS((PHITS_outSphere,PHITS_outBox))
       PHITSfile.writePHITS(phitsFilename)   
