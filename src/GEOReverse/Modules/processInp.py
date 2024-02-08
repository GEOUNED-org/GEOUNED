import configparser
from GEOReverse.Modules.options import Options


def setOptions(optionFile):

   setting = {'fileIn'   : ''               ,
              'fileOut'  : ''               ,
              'outBox'   : tuple()          ,
              'inFormat' : 'mcnp'           ,

              'UStart'   : 0                ,
              'levelMax' : 'all'            ,
              'cell'     : ('all',)         ,
              'mat'      : ('exclude',(0,)) }


   sectionNames =('Setting', 'Levels', 'Cells', 'Materials','Options')

   config = configparser.ConfigParser()
   config.optionxform = str
   config.read(optionFile)


   fileData = False
   for section in config.sections():
       if section not in sectionNames:
            print('{} bad section name'.format(section))
            exit()

       if section == 'Setting' :
           fileIn,fileOut,outBox,inFormat = getSetting(config)

           i = 0
           if fileIn  is not None : 
              setting['fileIn']  = fileIn
              i += 1

           if fileOut is not None : 
              setting['fileOut'] = fileOut
              i += 1

           if outBox  is not None : 
              setting['outBox']  = outBox
              i += 1

           if inFormat  is not None : 
              setting['inFormat']  = inFormat

           if i==3 : fileData = True


       elif section == 'Levels' :
           UStart,levMax = getLevel(config)
           if UStart is not None: setting['UStart'] = UStart
           if levMax is not None: setting['levelMax'] = levMax


       elif section in ('Cells','Materials'): 
           CMRange = getRange(section,config)
           if CMRange is None : continue
           if section == 'Cells':
              setting['cell'] = CMRange
           else:
              setting['mat'] = CMRange

       elif section == 'Options' :
           setSecOptions(config)

   if not fileData:
       print('missing input data in [Setting] section')
       exit()


   return setting


def getSetting(config):

   fileIn, fileOut, outBox,inFormat = None, None, None, None
   for key in config['Setting'].keys() :

      if  key in ['inputFile','mcnpFile']:
         fileIn  =  config.get('Setting',key)
      elif key == 'CADFile':
         fileOut =  config.get('Setting',key)
      elif key == 'outBox':
         outBox  = getBoxData(config.get('Setting',key))
      elif key == 'inFormat':
         inFormat  = config.get('Setting',key)
      else:
          print('{} bad keyword. Ignored'.format(key))       

   return fileIn, fileOut, outBox, inFormat 

      
def getLevel(config):

   UStart, levMax = None, None
   for key in config['Levels'].keys() :
      if  key == 'levelMax':
         levMax  =  config.get('Levels',key)
         if levMax.lower() == 'all' : levMax = 'all'
         elif levMax.isnumeric()    : levMax = int(UStart)
         else:
            print('bad value for levMax')
      elif key == 'UStart':
         UStart =  config.getint('Levels',key)
      else:
          print('{} bad keyword. Ignored'.format(key))       

   return UStart, levMax 


def getRange(section,config):
   rType   = None
   rawRange = None
   for key in config[section].keys() :
      if  key == 'rangeType':
         rType  =  config.get(section,key)
      elif key == 'range':
         rawRange  =  config.get(section,key)
      else:
          print('{} bad keyword. Ignored'.format(key))       

   
   if rType is None :
       return None
   elif rType.lower() == 'all':
       return ('all',)

   elif rType.lower() == 'exclude' or rType.lower() == 'include':
       ctRange = getRangeData(rawRange)
       if ctRange is None :
          print('bad range in section {}'.format(section))
          return None 
       else:
          return (rType.lower(),ctRange)
   else:
       return None     


def getBoxData(string):
   data = tuple(map(float,string.split()))
   if len(data) != 6 :
      print('bad Outbox value number') 
      exit()

   elif( (data[0] > data[1]) or
         (data[2] > data[3]) or
         (data[4] > data[5]) )   :
      print('bad Outbox boundaries') 
      exit()

   else:
      return (10*data[0],10*data[2],10*data[4],10*data[1],10*data[3],10*data[5])

 
def getRangeData(rawData):
   if rawData is None : return None

   intervals = rawData.split(',')
   rangeValues = set()

   for i in intervals:
      if i.strip().isnumeric() :
         rangeValues.add(int(i))
      else:
         bounds = i.split(':')
         if len(bounds) != 2:
            print('bad range definition')
            return None
         else:
            if bounds[0].strip().isnumeric() and bounds[1].strip().isnumeric():
               i1,i2 = map(int,bounds)
               if i1 > i2 :
                  print('bad range definition')
                  return None
               else:
                  rangeValues.update( range(i1,i2+1) ) 
            else:
               print('bad range definition')
               return None

   rangeValues = list(rangeValues)
   rangeValues.sort()
   return tuple(rangeValues)

def setSecOptions(config):

   for key in config['Options'].keys() :
      if  key == 'splitTolerance':
         tolerance  =  config.getfloat('Options',key)
         setattr(Options,key ,tolerance)
      else:
         print('{} bad keyword. Ignored'.format(key))       

   return  


def rangeGenerator(intervals):
   for i in intervals:
      if type(i) is tuple :
         for v in range(i[0],i[1]+1):
             yield v
      else:
         yield i




