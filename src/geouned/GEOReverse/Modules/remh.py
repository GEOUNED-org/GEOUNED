import math
import re

#########################################
# define patterns to be found in string #
#########################################


# used in card_split
celmat = re.compile(
    r"(?P<cel>^ *\d+) +(?P<scnd>(\d+|like))", re.I
)  # identify first two entries on cell card (cell name, material)
grlnumber = re.compile(
    r"[-+]?(\d+\.\d+|\.\d+|\d+\.?)(e[+-]\d+)?", re.I
)  # identify a general number form signed integer, float or exponential
param = re.compile(
    r"((^|\n {5})[\(\):\-\+\d+\.\# ]*)([\*a-z])", re.I
)  # identity begining of the paramter part of the cell card
likebut = re.compile(r"but", re.I)  # identify likebut  card
trans = re.compile(r"trcl|fill= *\d+[ c\$\n]*\(,", re.I)  # identify tranformed card
unsignedint = re.compile(r"\d+")

# user in get_stat function
reword = re.compile(r"(\d+|\(|\)|:|\#)")  # word identification in cell line
compcell = re.compile(r"\#\d")  # identify hashcell complement operator


# used in Complementary operator function
number = re.compile(r"(?P<number>[-+]?\d+)")  # signed (+-) (or not) numbers
# leftp=re.compile(r"^ *(?P<left>[-\d\(\#])",re.M)                               # identify first valid character
leftp = re.compile(r"(?P<left>[-+\d\(\#])", re.I)  # identify first valid character
rightp = re.compile(r"(?P<right>[ c\$\n]*$)", re.I)  # identify last valid character
# interblk=re.compile(r"(?P<previous>\d)(?P<next>(( +| *(\$)?\n(C\n)* *)[-+]?\d))") # two numbers separated by blank (or newline or comments)
# intercls=re.compile(r"(?P<previous>\))(?P<next>(( *| *(\$)?\n(C\n)* *)[-+]?\d))") # closed parenthesis followed by number
# interopn=re.compile(r"(?P<previous>\d)(?P<next>(( *| *(\$)?\n(C\n)* *)\())")   # number followed by opened parenthesis
# intercop=re.compile(r"(?P<previous>\))(?P<next>(( *| *(\$)?\n(C\n)* *)\())")   # closed parenthesis followed by opened parenthesis
interblk = re.compile(
    r"(?P<previous>\d)(?P<next>(( +| *((\n *)?\$|\nC)*\n *)[-+]?\d))"
)  # two numbers separated by blank (or newline or comments)
intercls = re.compile(
    r"(?P<previous>\))(?P<next>(( *| *((\n *)?\$|\nC)*\n *)[-+]?\d))"
)  # closed parenthesis followed by number
interopn = re.compile(r"(?P<previous>\d)(?P<next>(( *| *((\n *)?\$|\nC)*\n *)\())")  # number followed by opened parenthesis
intercop = re.compile(
    r"(?P<previous>\))(?P<next>(( *| *((\n *)?\$|\nC)*\n *)\())"
)  # closed parenthesis followed by opened parenthesis
colonamp = re.compile(r"[:&]")  # colon or amperserand

# used for remove redundant parenthesis function
mostinner = re.compile(r"\([^\(^\)]*\)")  # identify most inner parentheses
bracketsemi = re.compile(r"[\]\[;]")  # square bracket or semicolon
blnkline = re.compile(r"^ *\n", re.M)  # identify blank line
contline = re.compile(r"\n {0,4}(?P<start>[^c^ ])", re.I)  # identify character other than 'C' in fisrt 5 columns
comdollar = re.compile(r"\n(?P<blnk> *)\$")  # identify dollar on 'blank line'
startgeom = re.compile(r"(?P<previous>^ *)(?P<start>[\-\+\d])")  # identify beginning of the geomtric part
endgeom = re.compile(r"(?P<last>\d)(?P<next> *((\n *)?\$|\nc)?(\n *)?$)", re.I)  # identify end of the geomtric part
# endgeom=re.compile(r"(?P<last>\d)(?P<next> *(\$|\nc)?(\n *)?$)",re.I)                      # identify end of the geomtric part

# other
rehash = re.compile(r"# *(\d+|\()")  # find beginning of complementary operator (both cell and surf)
parent = re.compile(r"[\(|\)]")  # position of open and close parenthesis (get_hashcell)
gline = re.compile(
    r"(^ ?[\(\):\-\+\d+\.\# ]+|\n {5}[\(\):\-\+\d+\.\# ]+)", re.I
)  # valid geometric part of the line       (remove/restore_comments)
comments = re.compile(r"((\n *)?\$|\n *c)", re.I)  # begining of comment part               (remove/restore_comments)
# comments=re.compile(r"\$|\n *c",re.I)                               # begining of comment part               (remove/restore_comments)
celtrsf = re.compile(r"TRCL *= *", re.I)
celuniverse = re.compile(r"U *= *", re.I)
celfill = re.compile(r"FILL *= *", re.I)
trfnumber = re.compile(r"([-+]?(\d+\.\d+|\.\d+|\d+\.?)(e[+-]\d+)?)|\)", re.I)  # search for general number or close bracket )
likemat = re.compile(r"MAT *= *(?P<mat>\d+)", re.I)  # identify material value on like but card
dollar = re.compile(r"\$.*\n", re.I)


############################################################
# Auxiliary functions used in regular expresion functions  #
############################################################
def remove_dollar(string):
    m = dollar.search(string)
    while m:
        string = string[0 : m.start()] + " " + string[m.end() :]
        m = dollar.search(string)
    return string


def redundant(m, geom):
    """check if the inner parentheses are redundant"""
    term = m.group()

    # Find first valid character at the left of the  parenthese
    hashsmb = False
    leftOK = True
    left = m.start() - 1
    while left > -1:
        if geom[left] in ("\n", "C", "$", " "):
            left -= 1
        else:
            if geom[left] not in ("(", ":"):
                leftOK = False
            if geom[left] == "#":
                hashsmb = True
            break

    # if hash symbol found means parentheses delimits complementary
    # cell defined with surface. Theses parentheses are not redundants
    if hashsmb:
        return False

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


# function used in Regular expresion sub function
# function user in complementary function
# change the sign of the number
def chgsign(m):
    num = m.group(0)
    if num[0] == "-":
        return num[1:]
    if num[0] == "+":
        return "-" + num[1:]
    else:
        return "-" + num


# function used in Regular expersion sub function
# function user in complementary function
# change ':' in ')('  and
#        '&' in ':'
def repl_inter_union(m):
    if m.group(0) == ":":
        return ")("
    else:
        return ":"


# function used in Regular expersion sub function
# function user in remove_redundant function
# restore curve parentheses and colon characters
def reverse_repl(m):
    symb = m.group(0)
    if symb == "[":
        return "("
    elif symb == "]":
        return ")"
    else:
        return ":"


############################################################


def complementary(ccell, outter=True):
    """return the complementary cell"""
    wrkcell = Cline(ccell.str)
    if wrkcell.str[-1] == "\n":
        wrkcell.str = wrkcell.str[:-1]

    # simplify comment in geometry string
    wrkcell.remove_comments()

    # insert external parenthesis
    wrkcell.str = re.sub(leftp, r"(\g<left>", wrkcell.str, count=1)
    wrkcell.str = re.sub(rightp, r")\g<right>", wrkcell.str, count=1)

    # insert '&' as intersection operator
    wrkcell.str = re.sub(
        interblk, r"\g<previous>&\g<next>", wrkcell.str
    )  # change intersection separate by blank space ie: "number number"
    wrkcell.str = re.sub(interblk, r"\g<previous>&\g<next>", wrkcell.str)  # 2nd pass intersection blank space (require 2 pass)
    wrkcell.str = re.sub(
        intercls, r"\g<previous>&\g<next>", wrkcell.str
    )  # change intersection close parenthesis ie: ") number"
    wrkcell.str = re.sub(
        interopn, r"\g<previous>&\g<next>", wrkcell.str
    )  # change intersection open  parenthesis ie: "number ("
    wrkcell.str = re.sub(
        intercop, r"\g<previous>&\g<next>", wrkcell.str
    )  # change intersection close-open  parenthesis ie: ") ("

    # invert operators
    # substitute colon by ')(' and  '&' by colon
    wrkcell.str = re.sub(colonamp, repl_inter_union, wrkcell.str)
    wrkcell.remove_redundant(remove_com=False)

    # Change signs
    wrkcell.str = re.sub(number, chgsign, wrkcell.str)

    # insert external parenthesis
    if outter:
        wrkcell.str = re.sub(leftp, r"(\g<left>", wrkcell.str, count=1)
        wrkcell.str = re.sub(rightp, r")\g<right>", wrkcell.str, count=1)

    # restore original comments
    wrkcell.restore_comments()
    return wrkcell.str


############################################################
class Cline:
    def __init__(self, line):
        self.str = line

    def copy(self):
        return Cline(self.str)

    def get_surfaces_numbers(self):
        s = set(unsignedint.findall(self.str))
        return tuple(map(int, s))

    def remove_multispace(self):
        self.str.strip()
        self.str = re.sub(r" +", " ", self.str)

    def remove_cr(self):
        self.str = re.sub(r"\n", " ", self.str)
        return

    def remove_comments(self, full=False):
        """Remove the text of the comment. The symbol 'C' or '$' is
        kept in the line if 'full' option is False (default)"""
        celltab = re.split(gline, self.str)
        cont = True
        while cont:
            try:
                celltab.remove("")
            except:
                cont = False

        self.__comtab__ = []
        for i, s in enumerate(celltab):
            c = comments.match(s)
            if c:
                if not full:
                    self.__comtab__.append(s)
                    celltab[i] = c.group()
                else:
                    celltab[i] = ""

        self.str = "".join(celltab)
        return

    def restore_comments(self):
        """Restore the text of the comment."""
        celltab = re.split(gline, self.str)
        cont = True
        while cont:
            try:
                celltab.remove("")
            except:
                cont = False

        j = 0
        for i, s in enumerate(celltab):
            c = comments.match(s)
            if c:
                celltab[i] = self.__comtab__[j]
                j += 1

        self.str = "".join(celltab)
        return

    def outer_terms(self):
        cgeom = Cline(self.str)

        cgeom.remove_comments(full=True)
        cgeom.remove_redundant()
        cgeom.remove_cr()
        geom = cgeom.str

        # Loop until no redundant parentheses are found
        cont = True

        while cont:
            # Loop over most inner parentheses
            pos = 0
            cont = False
            while True:
                m = mostinner.search(geom, pos)
                if not m:
                    break
                cont = True
                if redundant(m, geom):
                    # remove redundant parentheses
                    geom = geom[: m.start()] + " " + geom[m.start() + 1 : m.end() - 1] + " " + geom[m.end() :]
                else:
                    # replace no redundant parentheses by 0 and : by ;
                    zeros = "0" * (m.end() - m.start())
                    geom = geom[: m.start()] + zeros + geom[m.end() :]

                pos = m.end()

        if ":" in geom:
            terms = []
            pos = 0
            while True:
                newpos = geom.find(":", pos)
                if newpos == -1:
                    terms.append(cgeom.str[pos:].strip())
                    break
                terms.append(cgeom.str[pos:newpos].strip())
                pos = newpos + 1
            return (terms, "OR")
        else:
            terms = []
            pos = 0
            while True:
                m = number.search(geom, pos)
                if not m:
                    break
                terms.append(cgeom.str[m.start() : m.end()])
                pos = m.end()
            return (terms, "AND")

    def remove_redundant(self, remove_com=True, remopt="nochg"):
        """return cell without redundant parenthesis"""

        # simplify comment in geometry string
        if remove_com:
            self.remove_comments()
        geom = self.str

        if remopt == "nochg" and geom.find(")") == -1:
            self.removedp = None
            return

        porg = self.countP()
        # Loop until no redundant parentheses are found
        cont = True
        while cont:
            # Loop over most inner parentheses
            pos = 0
            cont = False
            while True:
                m = mostinner.search(geom, pos)
                if not m:
                    break
                cont = True
                if redundant(m, geom):
                    # remove redundant parentheses
                    geom = geom[: m.start()] + " " + geom[m.start() + 1 : m.end() - 1] + " " + geom[m.end() :]
                else:
                    # replace no redundant parentheses by [] and : by ;
                    term = geom[m.start() + 1 : m.end() - 1].replace(":", ";")
                    geom = geom[: m.start()] + "[" + term + "]" + geom[m.end() :]
                pos = m.end()

        # restore curved parenthesis and colon
        geom = re.sub(bracketsemi, reverse_repl, geom)

        # remove possible blank line
        geom = re.sub(blnkline, "", geom)

        # ensure 5 blanks continous line
        geom = re.sub(contline, r"\n     \g<start>", geom)

        if remopt != "all":
            # add parenthesis to set geom as MCNP complex cell
            if geom.find(":") == -1 and geom.find("#") == -1:
                geom = re.sub(startgeom, r"\g<previous>(\g<start>", geom)
                geom = re.sub(endgeom, r"\g<last>)\g<next>", geom)

        # restore original comments
        self.str = geom
        pmod = self.countP()
        if remove_com:
            self.restore_comments()

        # subtitute comment $ with  blank line
        self.str = re.sub(comdollar, r"\nC\g<blnk>", self.str)
        pdiff = [x - y for x, y in zip(pmod, porg)]
        self.removedp = pdiff
        return

    def get_hashcell(self, start=0):
        """get the complementary cell defined with surfaces combination"""
        count = 0
        for p in parent.finditer(self.str, start):
            if p.group() == "(":
                count += 1
            else:
                count -= 1
            if count == 0:
                end = p.end()
                cell = Cline(self.str[start + 1 : end])
                break
        return cell, end

    def replace(self, surf, new, pos=0):

        s1 = str(surf)
        s2 = str(new)
        m = unsignedint.search(self.str, pos)
        if not m:
            print("no number in line")
            return -1
        if surf == new:
            return m.start() + len(s2)
        elif m.group() == s1:
            self.str = self.str[0 : m.start()] + s2 + self.str[m.end() :]
            return m.start() + len(s2)
        else:
            print(f"number {surf} not found")
            return -1

    def countP(self):
        lp = self.str.count("(")
        rp = self.str.count(")")
        return (lp, rp)

    def SplitCell(self):
        self.remove_redundant(remopt="all")
        self.remove_comments()
        geom = self.str
        self.restore_comments()

        # insert '&' as intersection operator
        geom = re.sub(
            interblk, r"\g<previous>&\g<next>", geom
        )  # change intersection separate by blank space ie: "number number"
        geom = re.sub(interblk, r"\g<previous>&\g<next>", geom)  # 2nd pass intersection blank space (require 2 pass)
        geom = re.sub(intercls, r"\g<previous>&\g<next>", geom)  # change intersection close parenthesis ie: ") number"
        geom = re.sub(interopn, r"\g<previous>&\g<next>", geom)  # change intersection open  parenthesis ie: "number ("
        geom = re.sub(intercop, r"\g<previous>&\g<next>", geom)  # change intersection close-open  parenthesis ie: ") ("

        parts = []
        block = ""
        iopen = 0
        ctype = None

        for c in geom:
            if c == "\n":
                continue
            if c in [":", "&"]:
                if iopen == 0:
                    if ctype == None:
                        if c == "&":
                            ctype = "AND"
                        if c == ":":
                            ctype = "OR"
                    if ctype == "AND" and c == ":":
                        ctype = "OR"
                        break
            elif c == "(":
                iopen += 1
            elif c == ")":
                iopen -= 1

        for c in geom:
            if c == "\n":
                continue
            if c in [":", "&"]:
                if iopen == 0:
                    if ctype == "AND" and c == "&":
                        parts.append(block)
                        block = ""
                    if ctype == "OR" and c == ":":
                        parts.append(block)
                        block = ""
                else:
                    if c != "&":
                        block = block + c
            else:
                if c != "&":
                    block = block + c
                if c == "(":
                    iopen += 1
                elif c == ")":
                    iopen -= 1
        parts.append(block)
        return parts, ctype


############################################################
class CellCardString:

    def __init__(self, card):
        self.stat = {"word": None, "hashcell": None, "hashsurf": None, "hash": None}

        self.name = None
        self.MAT = None
        self.U = None
        self.FILL = None
        self.TR = None
        self.TRCL = None
        self.__card_split__(card)
        self.__get_data__()

    def __card_split__(self, cardin):
        """Split the card string in three parts :
           - headstr : string containing the cell name, mat number and density (if mat != 0) of the cell
           - geom    : Cline class containing the part of the geometric definition of the cell
           - param    : Cline class containing the cell parameters part
        hproc is true if the complementary operator of the cell can be substituted"""

        m = celmat.match(cardin)
        self.name = int(m.group("cel"))
        self.hproc = True
        self.likeCell = None
        if m.group("scnd").lower() == "like":
            self.headstr = cardin[: m.start("scnd")]
            s = likebut.search(cardin, m.end("scnd"))
            self.geom = Cline(cardin[m.start("scnd") : s.end()])
            self.parm = Cline(cardin[s.end() :])
            self.hproc = False
            mc = unsignedint.search(self.geom.str)
            self.likeCell = int(mc.group())
            ml = likemat.search(self.parm.str)
            if ml:
                self.MAT = int(ml.group("mat"))
        elif m.group("scnd") == "0":
            cstart = m.end("scnd")
            self.MAT = 0
        else:
            self.MAT = int(m.group("scnd"))
            p = grlnumber.search(cardin, m.end("scnd"))
            cstart = p.end()

        if self.hproc:
            self.headstr = cardin[:cstart]
            cellcard = Cline(cardin[cstart:])
            cellcard.remove_comments()
            m = param.search(cellcard.str)
            if m:
                linecut = cellcard.str[: m.end(1)].count("\n")
            else:
                linecut = cellcard.str.count("\n")

            cellcard.restore_comments()

            # look for the last line geometry string
            if linecut != 0:
                pos = 0
                c = 0
                while c != linecut:
                    pos = cellcard.str.find("\n", pos)
                    c += 1
                m = param.search(cellcard.str, pos)

            if m:
                start = m.end(1)
                self.geom = Cline(cellcard.str[:start])
                self.parm = Cline(cellcard.str[start:])

                # look for transformation in cell parameters
                self.parm.remove_comments()
                m = trans.search(self.parm.str)
                if m:
                    self.hproc = False
                self.parm.restore_comments()
            else:
                self.geom = Cline(cellcard.str)
                self.parm = Cline("")

        return

    def __get_data__(self):

        mt = celtrsf.search(self.parm.str)
        if mt:
            string = self.parm.str[mt.end() :]
            m = grlnumber.search(string)

            if self.parm.str[mt.start() - 1] == "*":
                angle = True
            else:
                angle = False

            if "(" in string[: m.end()]:
                self.TRCL = []
                pos = 0
                while True:
                    m = trfnumber.search(string, pos)
                    pos = m.end()
                    if not m or m.group() == ")":
                        break
                    self.TRCL.append(float(m.group()))
                if angle and len(self.TRCL) > 3:
                    self.TRCL[3:12] = list(map(math.radians, self.TRCL[3:12]))
                    self.TRCL[3:12] = list(map(math.cos, self.TRCL[3:12]))
            else:
                self.TRCL = int(m.group())

        # get universe number
        m = celuniverse.search(self.parm.str)
        if m:
            m = grlnumber.search(self.parm.str[m.end() :])
            self.U = int(m.group())

        # get Fill number
        mf = celfill.search(self.parm.str)
        if mf:
            fillstring = remove_dollar(self.parm.str[mf.end() :])
            m = grlnumber.search(fillstring)
            self.FILL = int(m.group())

            # look for transformation after fill card
            m = re.search(r" *\(", fillstring)
            if m:
                if self.parm.str[mf.start() - 1] == "*":
                    angle = True
                else:
                    angle = False

                fillstring = fillstring[m.end() :]
                m = re.search(r"\)", fillstring)
                string = fillstring[: m.start()]
                tr = string.split()
                if len(tr) == 1:
                    self.TR = int(tr[0])
                else:
                    self.TR = list(map(float, tr))
                    if angle and len(self.TR) > 3:
                        self.TR[3:12] = list(map(math.radians, self.TR[3:12]))
                        self.TR[3:12] = list(map(math.cos, self.TR[3:12]))

    def get_stat(self, remove_com=True):
        """Count and return the number of words and hashes on the line."""

        if remove_com:
            self.geom.remove_comments()

        words = len(reword.findall(self.geom.str))
        hashcell = len(compcell.findall(self.geom.str))
        hashtot = self.geom.str.count("#")

        self.stat = {
            "words": words,
            "hash": hashtot,
            "hascell": hashcell,
            "hashsur": hashtot - hashcell,
        }
        if remove_com:
            self.geom.restore_comments()
        return [words, hashtot, hashcell, hashtot - hashcell]

    def get_lines(self):
        """split string card in the format Cards.lines of Cards object"""
        card = self.headstr + self.geom.str + self.parm.str

        # remove blank line introduced during the process
        card = re.sub(blnkline, "", card)

        # subtitute comment $ with  blank line
        card = re.sub(comdollar, r"\nC\g<blnk>", card)

        # ensure 5 blanks continous line
        card = re.sub(contline, r"\n     \g<start>", card)

        if card[-1] == "\n":
            card = card[:-1]
        return list(map(lambda x: x + "\n", card.split("\n")))


def remove_hash(cards, cname, keepComments=True):
    def remove(card, cname, keepComments):
        """remove complementary operator and subtitute by complementary cell"""
        if "parser.Card" in str(type(card)):
            celline = "".join(card.lines)
            cardstr = CellCardString(celline)
        else:
            cardstr = card
        cardstr.get_stat()
        cardstr.geom.remove_comments(full=not keepComments)
        if (not cardstr.hproc) or (cardstr.stat["hash"] == 0):
            return cardstr.geom  # no complementary operator or cannot be # cannot be removed
        cell = Cline(cardstr.geom.str)
        # find all operators in the cell and
        # substitute all complementary operators
        # locate object in list to reverse iteration
        hashgroup = []
        start = 0

        lencel = len(cell.str)
        while True:
            ic = cell.str.lower().find("c", start)
            idol = cell.str.find("$", start)
            if idol < 0:
                idol = lencel
            if ic < 0:
                ic = lencel
            end = min(idol, ic)
            for m in rehash.finditer(cell.str, start, end):
                hashgroup.append(m)
            start = cell.str.find("\n", end)
            if end == lencel:
                break

        for m in reversed(hashgroup):
            start = m.start()
            if m.group(1) == "(":  # complementary cell defined as surface intersections
                hcell, end = cell.get_hashcell(start)
                cellmod = cell.str[0:start] + complementary(hcell) + cell.str[end:]
                cell = Cline(cellmod)
            else:
                hcname = int(m.group(1))  # complementary cell defined with other cell index
                newdef = remove(cards[hcname], hcname, keepComments)  # remove complementary operator in new cell if necessary
                end = m.end()
                cellmod = cell.str[0:start] + "      " + complementary(newdef) + " " + cell.str[end:]
                cell = Cline(cellmod)
        return cell

    newcell = remove(cards[cname], cname, keepComments)
    return Cline(newcell.str)
