import re
mostinner=re.compile(r"\([^\(^\)]*\)")                                      # identify most inner parentheses
number   =re.compile(r"(?P<value>[-+]?\d+)")
mix      =re.compile(r"(?P<value>([-+]?\d+|\[0+\]))")
TFX      =re.compile(r"(?P<value>[FTXo]+)")
PValue   =re.compile(r"P\d+")
NValue   =re.compile(r"N\d+")
conversion = {'T':True,'F':False,'X':None}

class BoolSequence:
    def __init__(self,definition=None,operator=None) :
        if definition :
            self.elements = []
            self.setDef(definition)
        else:    
            self.elements = []
            self.operator = operator
            self.level    = 0

    def __str__(self):
       out='{}['.format(self.operator)
       if type(self.elements) is bool : return ' True ' if self.elements else ' False ' 
       for e in self.elements:
          if type(e) is int or type(e) is bool or type(e) is str : 
             out += ' {} '.format(e)
          else :
             out += e.__str__()

       out += '] '
       return out

    def append(self,*seq):
       for s in seq:
         if type(s) is int :
             level = -1
             if s in self.elements :
                 continue
             elif -s in self.elements:
                self.level = -1
                if self.operator == 'AND' :
                   self.elements = False 
                else:
                   self.elements = True
                return
         elif type(s) is bool :
             if self.operator == 'AND' and     s  or\
                self.operator == 'OR'  and not s   : 
                continue
             else:
                self.elements = s
                self.level =  -1
                return
         else:
             level = s.level 
             if type(s.elements) is bool:
               if self.operator == 'AND' and not s.elements  or\
                  self.operator == 'OR'  and     s.elements     : 
                  self.level =  -1
                  self.elements = s.elements 
                  return
               else:
                  continue

         self.elements.append(s)
         self.level = max(self.level,level+1)

    def assign(self,Seq):
        self.operator = Seq.operator
        self.elements = Seq.elements
        self.level    = Seq.level 

    def copy(self):
        cp = BoolSequence()
        cp.operator= self.operator
        cp.level = self.level
        if type(self.elements) is bool:
           cp.elements = self.elements
        else:
           for e in self.elements:
              if type(e) is int :
                 cp.elements.append(e)
              else:
                 cp.elements.append(e.copy())
        return cp        
        
    def getComplementary(self):

       c = BoolSequence(operator=self.compOperator())
       c.level = self.level

       if self.level == 0:
          for e in self.elements:
             c.elements.append(-e) 
          return c
       else:
         self.groupSingle()
         for e in self.elements:
            c.elements.append( e.getComplementary() )
         return c        

                
    def compOperator(self):
       if self.operator == 'AND': 
          return 'OR'
       else:
          return 'AND'

    def simplify(self,CT):
        if self.level > 0 :
           for seq in self.elements :
              seq.simplify(CT)
           self.clean()
           self.joinOperators()
           self.levelUpdate()

        if type(self.elements) is not bool and (self.level !=0 or len(self.elements) > 1) :
           levIn = self.level
           self.simplifySequence(CT)
           self.joinOperators()
           self.levelUpdate()
 
           if self.level > levIn: 
              self.simplify(CT)

            

    def simplifySequence(self,CT):
       surfNames = self.getSurfacesNumbers()
       if not surfNames : return

       #print(CT)
       #print('start',self)
       newNames = surfNames
       for valname in surfNames:
          #print('factorise',valname)
          if valname in newNames: 
             self.factorize(valname,CT)
             if type(self.elements) is bool: return
             newNames = self.getSurfacesNumbers()
       self.joinOperators()
       self.levelUpdate()
       

    def simplify_old(self,CT,loop=0):
       surfNames = self.getSurfacesNumbers()
       if not surfNames : return
       #print(CT)
       #print('start',self)
       newNames = surfNames
       simplified = False
       for valname in surfNames:
          #print('factorise',valname)
          if valname in newNames: 
             chg = self.factorize(valname,CT)
             simplified = simplified or chg
             newNames = self.getSurfacesNumbers()
       self.joinOperators()
       
       if self.level == 0 or (self.elements) is bool : return

       if loop == 0:
          ANDSeq = BoolSequence(operator=self.operator)
          ORSeq  = BoolSequence(operator=self.operator)
          for e in self.elements:
             if e.operator == 'AND':
                ANDSeq.append(e) 
             else:
                ORSeq.append(e) 

          ANDSeq.simplify(CT,loop=1)
          ORSeq.simplify(CT,loop=1)
          newSeq = BoolSequence(operator=self.operator)
          if ANDSeq.elements: newSeq.append(ANDSeq) 
          if ORSeq.elements : newSeq.append(ORSeq) 
          newSeq.joinOperators()
          self.assign(newSeq)
       else:
          for e in reversed(self.elements):
              e.simplify(CT,loop=0)
              if type(e.elements) is bool:
                 if self.operator == 'AND' and e.elements == False or \
                    self.operator == 'OR' and e.elements == True  : 
                    self.elements = e.elements
                    self.level = 0
                    break
                 else:
                    self.elements.remove(e)

          if self.elements == []:
             self.level = 0
             if self.operator == 'AND' : 
                self.elements = True 
             else:
                self.elements = False

     # check if level 0 sequence have oposite value a & -a = 0  , a|-a = 1
     # return the value of the sequence None(unknown), True, False
    def check(self):
       if self.level == 0:
          if type(self.elements) is bool : return self.elements

          signedSurf = set(self.elements)
          surfname  = self.getSurfacesNumbers()
          if len(signedSurf) == len(surfname) :  return None  # means same surface has not positive and negative value
          elif self.operator == 'AND' :          return False
          else:                                  return True
       else:
           if type(self.elements) is bool : return self.elements

           self.groupSingle()
           noneVal = False
           for e in self.elements:
              res = e.check()
              if res is None : noneVal = True
              elif self.operator == 'AND' and res is False : return False
              elif self.operator == 'OR'  and res is True  : return True

           if    noneVal :                return None
           elif  self.operator == 'AND' : return True
           else:                          return False


    def substitute(self,var,val):
       if   type(self.elements) is bool: return
       name = abs(var)
       ic = len(self.elements)
       for e in reversed(self.elements):
          ic -= 1
          if type(e) is int:
             if abs(e) == name :
                if type(val) is int : 
                   if name == e : self.elements[ic] = val
                   else         : self.elements[ic] = -val

                else :
                   if name == e : boolValue = val
                   else         : boolValue = not val

                   if self.operator == 'AND' and not boolValue :
                        self.elements = False
                        self.level    = 0
                        return
                   elif self.operator == 'OR' and boolValue :
                        self.elements = True
                        self.level    = 0
                        return
                   else:
                        self.elements.remove(e)
          else:
             e.substitute(var,val) 
 
       self.clean()
       self.levelUpdate()



    # remove sequence whom elements are boolean values instead of list
    def clean(self):   
        if type(self.elements) is bool : return self.elements
        for e in reversed(self.elements) :
           if type(e) is int : continue
           eVal = e.clean()
           if type(eVal) is bool:
              if eVal and self.operator == 'OR'  : 
                 self.elements = True
                 self.level = 0
                 return True
              elif not eVal and self.operator == 'AND' : 
                 self.elements = False
                 self.level = 0
                 return False
              self.elements.remove(e)

        if self.elements == [] :
           if self.operator == 'OR' : self.elements = False
           else                     : self.elements = True
           self.level = 0
           return self.elements
        else:
           return None 


    # join redundant operators in sequence
    def joinOperators(self):
        self.levelUpdate()
        if self.level == 0 : return
        if type(self.elements) is bool: return
        self.groupSingle()
        ANDop = []
        ORop  = []
     
        for e in self.elements:
           if e.operator == 'AND': ANDop.append(e) 
           else                  : ORop.append(e)

        if   len(ANDop) > 1  and self.operator == 'AND':
           newSeq = BoolSequence(operator='AND')
           for s in ANDop :
             newSeq.elements.extend(s.elements)
             self.elements.remove(s)
           newSeq.levelUpdate()
           self.append(newSeq)

            
        elif len(ORop)  > 1  and self.operator == 'OR':
           newSeq = BoolSequence(operator='OR')
           for s in ORop :
             newSeq.elements.extend(s.elements)
             self.elements.remove(s)
           newSeq.levelUpdate()
           self.append(newSeq)

        if self.level > 0  and len(self.elements)==1 :
           self.operator = self.elements[0].operator
           self.elements = self.elements[0].elements
           self.level -= 1
           self.joinOperators()
       
        if self.level == 0 : return
        for e in self.elements:
           e.joinOperators()



    def factorize(self,valname,CT=None):

        if CT is None :
           trueSet  = {abs(valname) : True }
           falseSet = {abs(valname) : False }
        else:
           trueSet,falseSet =  CT.getConstraintSet(valname)

        if trueSet is None:                    # valname cannot take True value 
             self.substitute(valname,False)
             return True

        if falseSet is None:                  # valname cannot take false value
             self.substitute(valname,True)
             return True

        funcVal = self.evaluate(trueSet)
        if funcVal == False :   
            newSeq = BoolSequence(operator='AND')
            #self.substitute(valname,False)
            for name,value in falseSet.items():
               self.substitute(name,value)
            newSeq.append(-valname,self.copy())
            self.assign(newSeq)
            return True

        elif funcVal == True:
            newSeq = BoolSequence(operator='OR')
            #self.substitute(valname,False)
            for name,value in falseSet.items():
               self.substitute(name,value)
            newSeq.append(valname,self.copy())
            self.assign(newSeq)
            return True


        funcVal = self.evaluate(falseSet)
        if funcVal == False :   
            newSeq = BoolSequence(operator='AND')
            #self.substitute(valname,True)
            for name,value in trueSet.items():
               self.substitute(name,value)
            newSeq.append(valname,self.copy())
            self.assign(newSeq)
            return True

        elif funcVal == True:
            newSeq = BoolSequence(operator='OR')
            #self.substitute(valname,True)
            for name,value in trueSet.items():
               self.substitute(name,value)
            newSeq.append(-valname,self.copy())
            self.assign(newSeq)
            return True

        return False
  

    def evaluate_newbad(self,valueSet,CT=None):

        if type(self.elements) is bool : return self.elements
        self.groupSingle()
        newSeq = self.copy()
        for name,value in valueSet.items():
            newSeq.substitute(name,value)
        
        if type(newSeq.elements) is bool :
            return newSeq.elements
        else :
            surfNames = tuple(newSeq.getSurfacesNumbers())        
            valname = surfNames[0]
            if CT is None :
                trueSet  = {abs(valname) : True }
                falseSet = {abs(valname) : False }
            else:
                trueSet,falseSet =  CT.getConstraintSet(valname)
                 
            if trueSet is None:
                trueVal = True
            else:       
                trueVal  = newSeq.evaluate(trueSet)
                if trueVal is None : return None
               
            if falseSet is None:
                falseVal = False
            else:       
                falseVal  = newSeq.evaluate(falseSet)
                if falseVal is None : return None
                    
            if trueVal != falseVal :
                return None
            else:
                return trueVal 


    def evaluate(self,valueSet):
        if type(self.elements) is bool : return self.elements
        self.groupSingle()
        op = self.operator
        noneVal = False
        for e in self.elements:
            if type(e) is int:
               key = abs(e)
               if key in valueSet.keys():
                    val = valueSet[key]
               else:
                    val = None
                    
               if val is not None : 
                  if e < 0 : val = not val
               else:
                  noneVal = True
            else:
               val = e.evaluate(valueSet)

            if val is None : noneVal = True
            elif op == 'AND' and val is False : return False
            elif op == 'OR'  and val is True  : return True

        if noneVal : 
           return None      
        if op == 'AND' :
           return True
        else:     
           return False
        
    def setDef(self,expression):
       terms,operator = outterTerms(expression)
       self.operator = operator
       self.level = 0
       lev0Seq    = set()
       lev0SeqAbs = set() 
       for t in terms :
         if isInteger(t) :
            val = int(t.strip('(').strip(')'))
            lev0Seq.add( val )
            lev0SeqAbs.add(abs(val))
            #self.elements.append(int(t.strip('(').strip(')')))
         else:
            x = BoolSequence(t) 
            self.level = max(x.level+1,self.level)
            self.append(x)
       
       # check if in integer sequence there is surface sequence s -s 
       if len(lev0Seq) != len(lev0SeqAbs) :
          if self.operator == 'AND' :
             self.elements = False
          else:
             self.elements = True
          self.level = -1
       else:
           self.append(*lev0Seq)
          
       self.groupSingle()     

    def groupSingle(self):
       if self.level == 0 : return
       if type(self.elements) is bool : return
       group = []
       for e in reversed(self.elements):
         if type(e) is int : 
            group.append(e)
            self.elements.remove(e)
         elif e.level==0 and len(e.elements) == 1 :
            group.append(e.elements[0])
            self.elements.remove(e)
            

       if not group : return
       seq = BoolSequence()
       seq.elements.extend(group)
       seq.operator = self.operator
       seq.level = 0
       self.elements.insert(0,seq)


    def getSurfacesNumbers(self):
        if type(self.elements) is bool : return tuple()
        surf = set()
        for e in self.elements:
           if type(e) is int :
                surf.add(abs(e))
           else:
                surf.update(e.getSurfacesNumbers())
        return surf             

    def levelUpdate(self):
       if type(self.elements) is bool :
          self.level = 0
          return
     
       self.level = 0
       for e in self.elements:
          if type(e) is  int : continue
          e.levelUpdate()
          self.level = max(e.level+1,self.level)

def insertInSequence(Seq,trgt,nsrf,operator):

    if operator == 'OR' :
        newSeq = BoolSequence(f'{trgt}:{nsrf}')
    else:
        newSeq = BoolSequence(f'{trgt} {nsrf}')

    substituteIntegerElement(Seq,trgt,newSeq)        
    Seq.joinOperators()


def substituteIntegerElement(Seq,target,newElement):
    for i,e in enumerate(Seq.elements):
       if type(e) is int:
          if e  == target :
               Seq.elements[i] = newElement
               Seq.level = max(Seq.level,1)
       else:
          substituteIntegerElement(e,target,newElement)
          

def outterTerms(expression,value='number'):
      if value == 'number' :
          #reValue = number
          reValue = mix
          nullVal = '0'
      else:
          reValue = TFX
          nullVal = 'o'
          
      expr = expression
      

      # Loop until no redundant parentheses are found
      cont = True
      
      while cont:
        # Loop over most inner parentheses
        pos = 0
        cont = False
        while True :
          m = mostinner.search(expr,pos)
          if not m : break
          cont = True
          if redundant(m,expr):
             # remove redundant parentheses
             expr = expr[:m.start()]+ ' ' + expr[m.start()+1:m.end()-1]+ ' ' + expr[m.end():]
          else:
             # replace no redundant parentheses by 0 and : by ;
             zeros = '[' + nullVal* (m.end()-m.start()-2) + ']' 
             expr = expr[:m.start()] + zeros + expr[m.end():]

          pos = m.end()

      if ':' in expr :
          terms = []
          pos = 0
          while True :
              newpos = expr.find(':',pos)
              if newpos == -1 :
                  terms.append(expression[pos:].strip())
                  break
              terms.append(expression[pos:newpos].strip())
              pos = newpos + 1                       
          return (terms,'OR')
      else:
          terms = []
          pos = 0
          while True:
              m = reValue.search(expr,pos)
              if not m : break
              terms.append(expression[m.start():m.end()])
              pos = m.end()                          
          return (terms,'AND')
        

def redundant(m,geom):
   """ check if the inner parentheses are redundant """
   term = m.group()

   # Find first valid character at the left of the  parenthese
   leftOK= True
   left = m.start()-1
   while left > -1:
       if geom[left] in ('\n','C','$',' '):
          left -= 1
       else:
          if geom[left] not in ('(',':') : leftOK  = False
          break

  # check if no ':' (or) are inside the parenthese
  # if not, parentheses are redundants
   if (term.find(':') == -1) : return True

  # Find first valid character at the right of the  parenthese
   rightOK= True
   right = m.end()
   while right < len(geom)  :
       if geom[right] in ('\n','C','$',' '):
          right += 1
       else:
          if geom[right] not in (')',':') : rightOK  = False
          break

  # if parentheses are like:
  # {( or : } ( ....... ) {) or :}
  # parentheses are redundants

   if leftOK and rightOK :
       return True
   else:
       return False

def isInteger(x):
    try :
      int(x.strip('(').strip(')'))
      return True
    except:
      return False  
