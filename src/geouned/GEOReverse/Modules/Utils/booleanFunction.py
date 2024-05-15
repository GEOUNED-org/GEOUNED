import re

mostinner = re.compile(r"\([^\(^\)]*\)")  # identify most inner parentheses
number = re.compile(r"(?P<value>[-+]?\d+)")
mix = re.compile(r"(?P<value>([-+]?\d+|\[0+\]))")
TFX = re.compile(r"(?P<value>[FTXo]+)")
PValue = re.compile(r"P\d+")
NValue = re.compile(r"N\d+")
conversion = {"T": True, "F": False, "X": None}


class BoolSequence:
    def __init__(self, definition=None, operator=None):
        if definition:
            self.elements = []
            self.set_def(definition)
        else:
            self.elements = []
            self.operator = operator
            self.level = 0

    def __str__(self):
        out = f"{self.operator}["
        if type(self.elements) is bool:
            return " True " if self.elements else " False "
        for e in self.elements:
            if type(e) is int or type(e) is bool or type(e) is str:
                out += f" {e} "
            else:
                out += e.__str__()

        out += "] "
        return out

    def append(self, *seq):
        for s in seq:
            if type(s) is int:
                level = -1
            else:
                level = s.level
                if type(s.elements) is bool:
                    if self.operator == "AND" and s.elements is False or self.operator == "OR" and s.elements is True:
                        self.level = 0
                        self.elements = s.elements
                        return
                    else:
                        continue

            self.elements.append(s)
            self.level = max(self.level, level + 1)

    def assign(self, Seq):
        self.operator = Seq.operator
        self.elements = Seq.elements
        self.level = Seq.level

    def copy(self):
        cp = BoolSequence()
        cp.operator = self.operator
        cp.level = self.level
        if type(self.elements) is bool:
            cp.elements = self.elements
        else:
            for e in self.elements:
                if type(e) is int:
                    cp.elements.append(e)
                else:
                    cp.elements.append(e.copy())
        return cp

    def get_complementary(self):

        c = BoolSequence(operator=self.comp_operator())
        c.level = self.level

        if self.level == 0:
            for e in self.elements:
                c.elements.append(-e)
            return c
        else:
            self.group_single()
            for e in self.elements:
                c.elements.append(e.get_complementary())
            return c

    def comp_operator(self):
        if self.operator == "AND":
            return "OR"
        else:
            return "AND"

    def simplify(self, CT, loop=0):
        surfNames = self.get_surfaces_numbers()
        if not surfNames:
            return
        #   print(CT)
        #    print('start',self)
        newNames = surfNames
        simplified = False
        for valname in surfNames:
            #      print('factorise',valname)
            if valname in newNames:
                chg = self.factorize(valname, CT)
                simplified = simplified or chg
                newNames = self.get_surfaces_numbers()
        self.join_operators()

        if self.level == 0 or (self.elements) is bool:
            return

        # new temp para ver como va
        # for e in self.elements:
        #   if type(e) is not int:
        #     e.simplify(CT)
        # return
        # end new

        if loop == 0:
            ANDSeq = BoolSequence(operator=self.operator)
            ORSeq = BoolSequence(operator=self.operator)
            for e in self.elements:
                if e.operator == "AND":
                    ANDSeq.append(e)
                else:
                    ORSeq.append(e)

            ANDSeq.simplify(CT, loop=1)
            ORSeq.simplify(CT, loop=1)
            newSeq = BoolSequence(operator=self.operator)
            if ANDSeq.elements:
                newSeq.append(ANDSeq)
            if ORSeq.elements:
                newSeq.append(ORSeq)
            newSeq.join_operators()
            self.assign(newSeq)
        else:
            for e in reversed(self.elements):
                e.simplify(CT, loop=0)
                if type(e.elements) is bool:
                    if self.operator == "AND" and e.elements is False or self.operator == "OR" and e.elements is True:
                        self.elements = e.elements
                        self.level = 0
                        break
                    else:
                        self.elements.remove(e)

            if self.elements == []:
                self.level = 0
                if self.operator == "AND":
                    self.elements = True
                else:
                    self.elements = False

    # check if level 0 sequence have oposite value a & -a = 0  , a|-a = 1
    # return the value of the sequence None(unknown), True, False
    def check(self):
        if self.level == 0:
            if type(self.elements) is bool:
                return self.elements

            signedSurf = set(self.elements)
            surfname = self.get_surfaces_numbers()
            if len(signedSurf) == len(surfname):
                return None  # means same surface has not positive and negative value
            elif self.operator == "AND":
                return False
            else:
                return True
        else:
            if type(self.elements) is bool:
                return self.elements

            self.group_single()
            noneVal = False
            for e in self.elements:
                res = e.check()
                if res is None:
                    noneVal = True
                elif self.operator == "AND" and res is False:
                    return False
                elif self.operator == "OR" and res is True:
                    return True

            if noneVal:
                return None
            elif self.operator == "AND":
                return True
            else:
                return False

    def substitute(self, var, val):
        if type(self.elements) is bool:
            return
        name = abs(var)
        ic = len(self.elements)
        for e in reversed(self.elements):
            ic -= 1
            if type(e) is int:
                if abs(e) == name:
                    if type(val) is int:
                        if name == e:
                            self.elements[ic] = val
                        else:
                            self.elements[ic] = -val

                    else:
                        if name == e:
                            boolValue = val
                        else:
                            boolValue = not val

                        if self.operator == "AND" and not boolValue:
                            self.elements = False
                            self.level = 0
                            return
                        elif self.operator == "OR" and boolValue:
                            self.elements = True
                            self.level = 0
                            return
                        else:
                            self.elements.remove(e)
            else:
                e.substitute(var, val)

        self.clean()
        self.level_update()

    def removeSurf(self, name):
        if type(self.elements) is bool:
            return
        ic = len(self.elements)
        for e in reversed(self.elements):
            ic -= 1
            if type(e) is int:
                if e == name:
                    if self.operator == "AND":
                        self.elements.remove(e)
                    elif self.operator == "OR":
                        self.elements = True
                        self.level = 0
                        return
            else:
                e.removeSurf(name)

        self.clean()
        self.level_update()

    # remove sequence whom elements are boolean values instead of list
    def clean(self):
        if type(self.elements) is bool:
            return self.elements
        for e in reversed(self.elements):
            if type(e) is int:
                continue
            eVal = e.clean()
            if type(eVal) is bool:
                if eVal and self.operator == "OR":
                    self.elements = True
                    self.level = 0
                    return True
                elif not eVal and self.operator == "AND":
                    self.elements = False
                    self.level = 0
                    return False
                self.elements.remove(e)

        if self.elements == []:
            if self.operator == "OR":
                self.elements = False
            else:
                self.elements = True
            self.level = 0
            return self.elements
        else:
            return None

    # join redundant operators in sequence
    def join_operators(self):
        if self.level == 0:
            return
        if type(self.elements) is bool:
            return
        self.group_single()
        ANDop = []
        ORop = []

        for e in self.elements:
            if e.operator == "AND":
                ANDop.append(e)
            else:
                ORop.append(e)

        if len(ANDop) > 1 and self.operator == "AND":
            newSeq = BoolSequence(operator="AND")
            for s in ANDop:
                newSeq.elements.extend(s.elements)
                self.elements.remove(s)
            newSeq.level_update()
            self.append(newSeq)

        elif len(ORop) > 1 and self.operator == "OR":
            newSeq = BoolSequence(operator="OR")
            for s in ORop:
                newSeq.elements.extend(s.elements)
                self.elements.remove(s)
            newSeq.level_update()
            self.append(newSeq)

        if self.level > 0 and len(self.elements) == 1:
            self.operator = self.elements[0].operator
            self.elements = self.elements[0].elements
            self.level -= 1
            self.join_operators()

        if self.level == 0:
            return
        for e in self.elements:
            e.join_operators()

    def factorize(self, valname, CT=None):

        if CT is None:
            trueSet = {abs(valname): True}
            falseSet = {abs(valname): False}
        else:
            trueSet, falseSet = CT.get_constraint_set(valname)

        if trueSet is None:
            self.substitute(valname, False)
            return True

        if falseSet is None:
            self.substitute(valname, True)
            return True

        funcVal = self.evaluate(trueSet, CT)
        if funcVal is False:
            newSeq = BoolSequence(operator="AND")
            # self.substitute(valname,False)
            for name, value in falseSet.items():
                self.substitute(name, value)
            newSeq.append(-valname, self.copy())
            self.assign(newSeq)
            return True
        elif funcVal is True:
            newSeq = BoolSequence(operator="OR")
            # self.substitute(valname,False)
            for name, value in falseSet.items():
                self.substitute(name, value)
            newSeq.append(valname, self.copy())
            self.assign(newSeq)
            return True

        funcVal = self.evaluate(falseSet, CT)
        if funcVal is False:
            newSeq = BoolSequence(operator="AND")
            # self.substitute(valname,True)
            for name, value in trueSet.items():
                self.substitute(name, value)
            newSeq.append(valname, self.copy())
            self.assign(newSeq)
            return True
        elif funcVal is True:
            newSeq = BoolSequence(operator="OR")
            # self.substitute(valname,True)
            for name, value in trueSet.items():
                self.substitute(name, value)
            newSeq.append(-valname, self.copy())
            self.assign(newSeq)
            return True

        return False

    def evaluate_newbad(self, valueSet, CT=None):

        if type(self.elements) is bool:
            return self.elements
        self.group_single()
        newSeq = self.copy()
        for name, value in valueSet.items():
            newSeq.substitute(name, value)

        if type(newSeq.elements) is bool:
            return newSeq.elements
        else:
            surfNames = newSeq.get_surfaces_numbers()
            valname = surfNames[0]
            if CT is None:
                trueSet = {abs(valname): True}
                falseSet = {abs(valname): False}
            else:
                trueSet, falseSet = CT.get_constraint_set(valname)

            if trueSet is None:
                trueVal = True
            else:
                trueVal = newSeq.evaluate(trueSet, CT)
                if trueVal is None:
                    return None

            if falseSet is None:
                falseVal = False
            else:
                falseVal = newSeq.evaluate(falseSet, CT)
                if falseVal is None:
                    return None

            if trueVal != falseVal:
                return None
            else:
                return trueVal

    def evaluate(self, valueSet, CT=None):
        if type(self.elements) is bool:
            return self.elements
        self.group_single()
        op = self.operator
        noneVal = False
        for e in self.elements:
            if type(e) is int:
                key = abs(e)
                if key in valueSet.keys():
                    val = valueSet[key]
                else:
                    val = None

                if val is not None:
                    if e < 0:
                        val = not val
                else:
                    noneVal = True
            else:
                val = e.evaluate(valueSet)

            if val is None:
                noneVal = True
            elif op == "AND" and val is False:
                return False
            elif op == "OR" and val is True:
                return True

        if noneVal:
            return None
        if op == "AND":
            return True
        else:
            return False

    def set_def(self, expression):
        terms, operator = outer_terms(expression)
        self.operator = operator
        self.level = 0
        for t in terms:
            if is_integer(t):
                self.elements.append(int(t.strip("(").strip(")")))
            else:
                x = BoolSequence(t)
                self.level = max(x.level + 1, self.level)
                self.elements.append(x)
        self.group_single()

    def group_single(self):
        if self.level == 0:
            return
        if type(self.elements) is bool:
            return
        group = []
        for e in reversed(self.elements):
            if type(e) is int:
                group.append(e)
                self.elements.remove(e)
            elif e.level == 0 and len(e.elements) == 1:
                group.append(e.elements[0])
                self.elements.remove(e)

        if not group:
            return
        seq = BoolSequence()
        seq.elements.extend(group)
        seq.operator = self.operator
        seq.level = 0
        self.elements.insert(0, seq)

    def sameOperator(self):
        if self.level == 0:
            return
        if type(self.elements) is bool:
            return

        changed = []
        new = []
        for i, e in enumerate(self.elements):
            if e.operator == self.operator:
                changed.append(i)
                for s in e.elements:
                    if type(s) is int:
                        new.append(BoolSequence(str(s)))
                    else:
                        new.append(s)

        for i in reversed(changed):
            del self.elements[i]
        self.append(*new)
        self.level_update()

        return

    def get_surfaces_numbers(self):
        if type(self.elements) is bool:
            return tuple()
        surf = set()
        for e in self.elements:
            if type(e) is int:
                surf.add(abs(e))
            else:
                surf.update(e.get_surfaces_numbers())
        return tuple(surf)

    def level_update(self):
        if type(self.elements) is bool:
            self.level = 0
            return

        newlev = 0
        for e in self.elements:
            if type(e) is int:
                lev = 0
            else:
                lev = e.level + 1
            newlev = max(lev, newlev)
        self.level = newlev


def outer_terms(expression, value="number"):
    if value == "number":
        # reValue = number
        reValue = mix
        nullVal = "0"
    else:
        reValue = TFX
        nullVal = "o"

    expr = expression

    # Loop until no redundant parentheses are found
    cont = True

    while cont:
        # Loop over most inner parentheses
        pos = 0
        cont = False
        while True:
            m = mostinner.search(expr, pos)
            if not m:
                break
            cont = True
            if redundant(m, expr):
                # remove redundant parentheses
                expr = expr[: m.start()] + " " + expr[m.start() + 1 : m.end() - 1] + " " + expr[m.end() :]
            else:
                # replace no redundant parentheses by 0 and : by ;
                zeros = "[" + nullVal * (m.end() - m.start() - 2) + "]"
                expr = expr[: m.start()] + zeros + expr[m.end() :]

            pos = m.end()

    if ":" in expr:
        terms = []
        pos = 0
        while True:
            newpos = expr.find(":", pos)
            if newpos == -1:
                terms.append(expression[pos:].strip())
                break
            terms.append(expression[pos:newpos].strip())
            pos = newpos + 1
        return (terms, "OR")
    else:
        terms = []
        pos = 0
        while True:
            m = reValue.search(expr, pos)
            if not m:
                break
            terms.append(expression[m.start() : m.end()])
            pos = m.end()
        return (terms, "AND")


def redundant(m, geom):
    """check if the inner parentheses are redundant"""
    term = m.group()

    # Find first valid character at the left of the  parenthese
    leftOK = True
    left = m.start() - 1
    while left > -1:
        if geom[left] in ("\n", "C", "$", " "):
            left -= 1
        else:
            if geom[left] not in ("(", ":"):
                leftOK = False
            break

    # check if no ':' (or) are inside the parenthese
    # if not, parentheses are redundants
    if term.find(":") == -1:
        return True

    # Find first valid character at the right of the  parenthese
    rightOK = True
    right = m.end()
    while right < len(geom):
        if geom[right] in ("\n", "C", "$", " "):
            right += 1
        else:
            if geom[right] not in (")", ":"):
                rightOK = False
            break

    # if parentheses are like:
    # {( or : } ( ....... ) {) or :}
    # parentheses are redundants

    if leftOK and rightOK:
        return True
    else:
        return False


def is_integer(x):
    try:
        int(x.strip("(").strip(")"))
        return True
    except:
        return False
