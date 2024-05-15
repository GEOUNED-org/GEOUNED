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
colonamp = re.compile(r"[:&]")  # colon or ampersand

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
lastCR = re.compile(r"\n *$")  # last newline "\n" on line string


# restore curve parentheses and colon characters
def reverse_repl(m):
    symb = m.group(0)
    if symb == "[":
        return "("
    elif symb == "]":
        return ")"
    else:
        return ":"


def countP(string):
    lp = string.count("(")
    rp = string.count(")")
    return (lp, rp)


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


def remove_redundant(geom):
    """return cell without redundant parenthesis"""

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

    # remove last newline at the end of the line string
    geom = re.sub(lastCR, "", geom)

    return geom
