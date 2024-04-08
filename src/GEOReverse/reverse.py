import sys
import FreeCAD, Part, Import

from GEOReverse.Modules.MCNPinput import MCNPinput
from GEOReverse.Modules.XMLinput import XMLinput
from GEOReverse.Modules.buildCAD import buildCAD,makeTree
from GEOReverse.Modules.Objects import CADCell
from GEOReverse.Modules.processInp import setOptions 
from GEOReverse.CodeVersion import *


def reverse(optFile='configRevese.ini') :

   printCodeVersion()

   setting = setOptions(optFile)

   geomfile = setting['fileIn']
   outname  = setting['fileOut']
   outBox   = setting['outBox']
   inFormat = setting['inFormat']

   CADselection = { 'Ustart'   : setting['UStart']  ,
                    'levelMax' : setting['levelMax'], 
                    'cell'     : setting['cell']    ,
                    'mat'      : setting['mat']     ,
                    'format'   : setting['inFormat']
                 }  

   UnivCell = CADCell()
   UnivCell.shape = UnivCell.makeBox(FreeCAD.BoundBox(*outBox))

   # get geometry definition from MCNP input
   if inFormat == 'mcnp':
      geom = MCNPinput(geomfile)
   elif inFormat == 'openMC_XML':
      geom = XMLinput(geomfile)
   else:
      msg = (
         f'input format type {inFormat} is not supported.'
         'Supported options are "mcnp" or "openMC_XML"'
      )
      raise ValueError(msg)

   CADCells,fails = buildCAD(UnivCell,geom,CADselection)

   if fails : print('failed in conversion',fails)

   CADdoc=FreeCAD.newDocument("WorkingDoc")

   makeTree(CADdoc,CADCells)
   Import.export(CADdoc.Objects[0:1],outname+'.stp' ) 
   CADdoc.saveAs(outname+'.FCStd')



def printCodeVersion():

   FreeCAD_Version = '{V[0]:}.{V[1]:}.{V[2]:}'.format(V=FreeCAD.Version())
   title = '''\
#########################################################################
#                                                                       # 
#      GEOReverse version {:<11}{}{:>26}
#      FreeCAD    version {:<11}{:>36}  
#                                                                       # 
#########################################################################'''.format(GEOReverse_Version,GEOReverse_ReleaseDate,'#',FreeCAD_Version,'#')
   print(title)


if __name__ == '__main__' :
  reverse()



