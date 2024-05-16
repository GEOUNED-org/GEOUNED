# -*- coding: utf-8 -*-

"""
Functions for parsing MCNP input files.
"""

import re
import warnings

from .PartialFormatter import PartialFormatter

# integer with one prefix character
re_int = re.compile(r"\D{0,1}\d+")

# interior of square brackets for tally in lattices
re_ind = re.compile(r"\[.+\]", flags=re.DOTALL)

# repitition syntax of MCNP input file
re_rpt = re.compile(r"\d+[ri]", flags=re.IGNORECASE)

# imp or tmp parameters in cell card
re_prm = re.compile(r"((imp:n|imp:p|tmp)\s+\S+)")
re_prm = re.compile(r"[it]mp:*[npe]*[=\s]+\S+", flags=re.IGNORECASE)
re_prm = re.compile(r"([it]mp:*[npe]*[=\s]+)(\S+)", flags=re.IGNORECASE)

# fill keyword
re_fll = re.compile(r"\*{0,1}fill[=\s]+", flags=re.IGNORECASE)  # TODO: this will also match fill===


# If type specifier not given, any data type can be formatted:
def fmt_gen(s):
    return "{" + f":<{len(s)}" + "}"


fmt_d = fmt_gen
fmt_g = fmt_gen
fmt_s = fmt_gen

partial_formmatter = PartialFormatter()


class CidClass(object):
    """
    There are two levels of card types. 1-st level is purely defined by card
    position in the input file.  There can be:
          * message card
          * title card
          * cell card
          * surface card
          * data card.

    Data cards can be of different type that is defined by the card's first
    entry. Therefore data cards can be characterized by type of the 2-nd level:
          * m cards,
          * f cards,
          * etc.

    This class is to describe the 1-st level of card types.
    """

    # no-information cards
    comment = -1
    blankline = -2
    # card types in order they appear in MCNP input
    message = 1
    title = 2
    cell = 3
    surface = 4
    data = 5

    @classmethod
    def get_name(cls, cid):
        """
        Return the name of the card type by its index.
        """
        for k, v in list(cls.__dict__.items()):
            if "__" not in k and v == cid:
                return k
        else:
            print("No attribute with name", cid)
            raise ValueError()


CID = CidClass()


class Card(object):
    """
    Representation of a card.
    """

    def __init__(self, lines, ctype, pos, debug=None):

        # Original lines, as read from the input file
        self.lines = lines

        # card type by its position in the input. See CID class.
        self.ctype = ctype

        # True if self.lines has changed after initialization
        # used in remove_hash function
        self.cstrg = False

        # data card type. Defined from the get_values() method.
        # Has sense only to data cards (see ctype). For other card types
        # is None.
        self.dtype = None

        # Input file line number, where the card was found.
        self.pos = pos

        # File-like object to write debug info (if not None)
        self.debug = debug

        # template string. Represents the general structure of the card. It is
        # a copy of lines, but meaningful parts are replaced by format
        # specifiers, {}
        self.template = ""

        # List of strings represenging meaningful parts of the card. The
        # original multi-line string card is obtained as
        # template.format(*input)
        self.input = []

        # Dictionary of parts that are removed from input before processing it.
        # For example, repitition syntax (e.g. 5r or 7i) is replaced with '!'
        # to prevent its modification.
        self.hidden = {}

        # List of (v, t) tuples, where v -- value and t -- its type.
        self.values = []

        # geometry prefix and suffix
        # self.geom_prefix = ''
        # self.geom_suffix = ''

        # some properties defined on demand
        # cell properties
        self.__u = -1  # -1 means undefined. None -- not specified in input
        self.__f = -1  # fill
        self.__m = -1  # material
        self.__d = ""  # density
        self.__i = -1  # importances
        self.__cr = -1  # set of reference cells.
        # surface properties
        self.__st = ""  # '' means undefined.

        # Split card to template and meaningful part is always needed. Other
        # operations are optional.
        self.get_input()
        return

    def _get_value_by_type(self, t):
        """
        Returns the first value of type t found in self.values.
        """
        vl, tl = zip(*self.values)
        try:
            i = tl.index(t)
        except ValueError:
            return None
        finally:
            return vl[i]

    def _set_value_by_type(self, t, v):
        """
        Sets the first value of type t to v in self.values.
        """
        vl, tl = zip(*self.values)
        i = tl.index(t)
        self.values[i] = (v, t)

    @property
    def geom_prefix(self):
        return self._get_value_by_type("#gpr")

    @geom_prefix.setter
    def geom_prefix(self, value):
        return self._set_value_by_type("#gpr", value)

    @property
    def geom_suffix(self):
        return self._get_value_by_type("#gsu")

    @geom_suffix.setter
    def geom_suffix(self, value):
        return self._set_value_by_type("#gsu", value)

    def print_debug(self, comment, key="tihv"):
        d = self.debug
        if d:
            print(
                f"Line {self.pos}, {CID.get_name(self.ctype)} card. {comment}",
                file=d,
            )
            if "t" in key:
                print("    template:", repr(self.template), file=d)
            if "i" in key:
                print("    input:   ", self.input, file=d)
            if "h" in key:
                print("    hidden:  ", self.hidden, file=d)
            if "v" in key:
                print("    values:  ", self.values, file=d)

    def get_input(self, check_bad_chars=False):
        """
        Recompute template, input and hidden attributes from lines
        """

        mline = "".join(self.lines)
        if check_bad_chars:
            bad_chars = "\t"
            for c in bad_chars:
                if c in mline:
                    if self.debug:
                        self.print_debug("get_input: bad char in input cards", "")
                    else:
                        raise ValueError("Bad character in input file. " + "Run with --debug option.")

        if self.ctype in (CID.comment, CID.blankline):
            # nothing to do for comments or blanklines:
            self.input = ""
            self.template = mline

        else:
            # TODO: protect { and } in comment parts of the card.
            tmpl = []  # part of template
            inpt = []  # input, meaningful parts of the card.
            if mline.split()[0][:2] == "fc":
                # this is tally comment. It always in one line and is not
                # delimited by & or $
                i = mline[:80]
                t = mline.replace(i, "{}", 1)
                inpt = [i]
                tmpl = [t]
            else:
                for l in self.lines:
                    if is_commented(l):
                        tmpl.append(l)
                    else:
                        # entries optionally delimited from comments by $ or &
                        # requires that delimiters prefixed with space
                        d = index_(l, "$&")
                        # d1 = l.find(' $')
                        # d2 = l.find(' &')
                        # if -1 < d1 and -1 < d2:
                        #     # both & and $ in the line. Use the most left
                        #     d = min(d1, d2)
                        # elif d1 == -1 and d2 == -1:
                        #     # no delimiters at all, whole line is meaningful
                        #     # except the new-line char
                        #     d = len(l) - 1
                        # else:
                        #     # only one of delimiters is given.
                        #     d = max(d1, d2)
                        i = l[:d]
                        t = l[d:]
                        inpt.append(i)
                        tmpl.append(fmt_s(i) + t)
            self.input = inpt
            self.template = "".join(tmpl)

            # TODO: dtype and name of the card can be defined already here.

        self.print_debug("get_input", "ti")
        return

    def _protect_nums(self):
        """
        In the meaningful part of the card replace numbers that do not
        represent cell, surface or a cell parameter with some unused char.
        """

        inpt = "\n".join(self.input)

        d = {}

        # in cell card:
        if self.ctype == CID.cell:  # and 'like' not in inpt:
            d["~"] = []  # float values in cells

            # Replace material density
            if "like" not in inpt:
                tokens = inpt.replace("=", " ").split()
                cell, mat, rho = tokens[:3]
                if int(mat) != 0:
                    for s in (cell, mat, rho):
                        inpt = inpt.replace(s, "~", 1)
                    inpt = inpt.replace("~", cell, 1)
                    inpt = inpt.replace("~", mat, 1)
                    d["~"].append(rho)

            # imp and tmp parameters:
            # print 're_prm: inp', repr(inpt)
            for s1, s2 in re_prm.findall(inpt):
                # print 're_prm: fnd', repr(s)
                d["~"].append(s2)
                inpt = inpt.replace(s1 + s2, s1 + "~", 1)

        # replace repitition syntax in junks:
        sbl = re_rpt.findall(inpt)
        if sbl:
            for s in sbl:
                inpt = inpt.replace(s, "!", 1)
            d["!"] = sbl

        if self.ctype == CID.data and inpt.lstrip().lower()[0] == "f" and inpt.lstrip()[1].isdigit():
            # this is tally card. Hide indexes in square brackets
            sbl = re_ind.findall(inpt)
            if sbl:
                for s in sbl:
                    inpt = inpt.replace(s, "|", 1)
                d["|"] = sbl

        self.input = inpt.split("\n")
        self.hidden = d

        self.print_debug("_protect_nums", "ih")
        return

    def get_values(self):
        """
        Replace integers in the meaningfull part with format specifiers, and
        populate the `values` attribute.
        """
        self._protect_nums()
        if self.ctype == CID.cell:
            inpt, vt = _split_cell(self.input, self)
            self.name = vt[0][0]
        elif self.ctype == CID.surface:
            inpt, vt, stype, scoef = _split_surface(self.input)
            self.TR = None
            if len(vt) > 1:
                if vt[1][1] == "tr":
                    self.TR = vt[1][0]
            self.stype = stype.lower()
            self.scoefs = scoef
            self.name = vt[0][0]
        elif self.ctype == CID.data:
            inpt, vt, dtype = _split_data(self.input)
            self.dtype = dtype
            if dtype == "TRn":
                # print(vt,inpt)
                unit, inpt, fvals = _parse_tr(inpt)
                self.unit = unit
                vt += fvals
            if self.dtype is not None:
                self.name = vt[0][0]
        else:
            inpt = self.input
            vt = []

        self.input = inpt
        self.values = vt

        self.print_debug("get_values", "iv")
        return

    def get_refcells(self):
        """
        Returns all cells used in definition of self.
        """
        if self.ctype != CID.cell:
            return None
        if self.__cr != -1:
            return self.__cr
        else:
            s = set()
            for v, t in self.values:
                if t == "cel":
                    s.add(v)
            self.__cr = s
            return self.__cr

    def get_geom(self):
        """
        Returns part of the cell card describing geometry, as a (multiline)
        string.
        """
        p, s = self.geom_prefix, self.geom_suffix
        self.geom_prefix = "§"
        self.geom_suffix = "§"
        geom = self.card().split("§")[1]
        self.geom_prefix = p
        self.geom_suffix = s
        return geom

    def get_u(self):
        """
        Returns universe, the cells belongs to.
        """
        if self.ctype != CID.cell:
            return None
        if self.__u != -1:
            return self.__u
        else:
            # get it only once:
            for v, t in self.values:
                if t == "u":
                    self.__u = v
                    break
            else:
                self.__u = None
            return self.__u

    def get_m(self):
        """
        For cell card return material number
        """
        if self.ctype != CID.cell:
            return None

        if self.__m != -1:
            return self.__m
        else:
            if "like" in "".join(self.input).lower():
                # material name should be given in another cell.
                pass
            for v, t in self.values:
                if t == "mat":
                    self.__m = v
                    break
            else:
                # raise ValueError("Cell does not have material specs")
                self.__m = -2
            return self.__m

    def get_d(self):
        """
        For cell card return density
        """
        if self.__d != "":
            return self.__d

        if self.get_m() == 0:
            self.__d = 0
            return self.__d
        elif self.get_m() == -2:
            # this is like-but cell
            self.__d = -100.0
            return self.__d
        else:
            # density entry is hidden in the input and available as the 1-st
            # entry in self.hidden dictionary.
            self.__d = float(self.hidden["~"][0])
            return self.__d

    def set_d(self, v):
        """
        Set density. Accespted are string represetntaion of a float.

        It is assumed that get_values() method is called before this.
        """
        if self.get_m() > 0:
            self.hidden["~"][0] = v
            self.__d = float(v)

    def get_f(self, newv=None, trsf=False):
        """
        Returns universe, the cell is filled with.
        """
        if self.ctype != CID.cell:
            return None

        if self.__f != -1 and newv is None:
            if trsf:
                params = []
                for v, t in self.values:
                    if t == "#tunit":
                        unit = v
                    elif t == "tr":
                        params = [v]
                        break
                    elif t == "#tparam":
                        params.append(float(v))
                return self.__f, (unit, params)
            else:
                return self.__f
        else:
            # get it only once:
            unit = ""
            params = []
            for i in range(len(self.values)):
                v, t = self.values[i]
                if t == "fill":
                    if newv is not None:
                        v = newv
                        self.values[i] = (v, t)
                    self.__f = v
                    if trsf:
                        unit = self.values[i - 1][0]
                        for v, t in self.values[i + 1 :]:
                            if t == "tr":
                                params = [v]
                                break
                            elif t == "#tparam":
                                params.append(float(v))
                    break
            else:
                self.__f = None

            if trsf:
                return self.__f, (unit, params)
            else:
                return self.__f

    def get_imp(self, vals={}):
        """
        Returns importances, if explicitly specified in the cell card.
        """
        if self.ctype != CID.cell:
            return None

        if self.__i != -1 and not vals:
            return self.__i
        else:
            res = {}
            inpt = " ".join(self.input).lower()
            for p in "npe":
                key = "imp:" + p

                s = inpt.split(key)
                if len(s) == 1:
                    # there is no key in the input line.
                    continue
                else:
                    n = s[0].count("~")
                    res[key] = float(self.hidden["~"][n])
                    if p in vals:
                        # change value only if necessary
                        if res[key] != vals[p]:
                            res[key] = vals[p]
                            self.hidden["~"][n] = str(vals[p])

                # for s in self.hidden.get('~', []):
                #     sl = s.lower()
                #     if key in sl:
                #         val = float(sl.replace(key, '').replace('=', ''))
                #         res[key] = val
            if not res:
                res["imp:n"] = 1
            self.__i = res
            return self.__i

    def remove_fill(self):
        """
        Removes the FILL= keyword of a cell card.

        This method must be called after get_values().
        """
        # Fill card is followed by one universe number and optionally by
        # transformation in parentheses. Optionally, the `fill` keyword can be
        # prefixed with asterisk

        # Two possibilites are here: (1) start from original line and clean out
        # everything related to FILL, or (2) modify already existing template,
        # values and input. In both cases, the new card should be checked for
        # empty lines.

        # Case (2): Modify existing input and values. All content after the
        # FILL keyword is parsed and thus is given in values, while input
        # provides placefor them. THus, simply replaceing values with spaces
        # will almost do the job. The remaining part -- the keyword itself that
        # is presented in input.

        # replace with spaces all FILL-related tokens
        vals = []  # new values list.
        oldv = self.values[:]
        state = "before"
        while oldv:
            v, t = oldv.pop(0)
            if state == "before" and t == "fill":
                v = " "
                state = "afterU"
            elif state == "afterU" and "(" in t:
                v = " "
                state = "after("
            elif state == "after(":
                v = " "
                if ")" in t:
                    state = "after"
            vals.append((v, t))
        self.values = vals

        # Remove FILL from the input
        for n, i in enumerate(self.input):
            if "fill" in i.lower():
                # This part of input contains the fill keyword. This keyword is
                # optionally prepended with an asterix and followed by a sign
                i = re_fll.sub(" ", i)
                self.input[n] = i
                break

        self.print_debug("remove_fill", "iv")
        return

    def card(self, wrap=False, comment=True):
        """
        Return multi-line string representing the card.
        """
        if self.input:
            # put values back to meaningful parts:
            inpt = "\n".join(self.input)
            inpt = inpt.format(*[t[0] for t in self.values])

            # put back hidden parts:
            for k, vl in list(self.hidden.items()):
                for v in vl:
                    inpt = inpt.replace(k, v, 1)

            inpt = inpt.split("\n")
            if not comment:
                return " ".join(inpt)

            if wrap:  # and self.ctype != CID.title:
                indent = " " * 5
                if self.ctype == CID.title:
                    indent = "c" + indent
                tparts = re.split(r"\{.*?\}", self.template)[1:]
                # print 'wrapped inp', repr(self.template)
                # print 'wrapped spl', repr(tparts)
                newt = [""]  # new template parts
                newi = []  # new input parts
                self.print_debug("card wrap=True", "")
                for i, t in zip(inpt, tparts):
                    self.print_debug("    " + repr(i) + repr(t), "")
                    il = []
                    tl = [t]

                    # while len(i.rstrip()) > 79:
                    while len(i.rstrip()) > 80:
                        # first try to shift to left
                        if i[:5] == " " * 5:
                            i = " " * 5 + i.lstrip()
                        if len(i.rstrip()) > 79:
                            # input i must be wrapped. Find proper place:
                            for dc in " :":
                                k = i.rstrip().rfind(dc, 0, 75)
                                if k > 6:
                                    il.append(i[:k])
                                    tl.append("\n")
                                    i = indent + i[k:]
                                    self.print_debug("card wrap=True" + repr(il[-1]) + repr(i), "")
                                    break
                            else:
                                # there is no proper place to wrap.
                                self.print_debug("Cannot wrap line " + repr(i), "")
                                warnings.warn(f"Cannot wrap card on line {self.pos}")
                                break
                        else:
                            # input i fits to one line. Do nothing.
                            pass

                    newt += tl
                    newi += il + [i]
                tmpl = "{}".join(newt)
                inpt = newi
            else:
                tmpl = self.template

            card = partial_formmatter.format(tmpl, *inpt)
            # card = tmpl.format(*inpt)
        else:
            card = self.template
        return card

    def remove_spaces(self):
        """
        Remove extra spaces from meaningful parts.
        """
        self.print_debug("before remove_spaces", "i")
        if self.ctype in (CID.cell, CID.surface, CID.data):
            inpt = []
            for i in self.input:
                indented = i[:5] == " " * 5
                # leave only one sep. space
                i = " ".join(i.split())
                i = i.strip()
                # spaces before/after some characters are not needed:
                for c in "):":
                    i = i.replace(" " + c, c)
                for c in "(:":
                    i = i.replace(c + " ", c)
                if indented:
                    i = " " * 5 + i
                inpt.append(i)
                self.print_debug(i, "")
            self.input = inpt
            self.print_debug("after remove_spaces", "i")
        return

    def apply_map(self, f):
        """
        Replace Ni in self.values by Mi = f(Ni, Ti).
        """
        self.print_debug("before apply_map", "vi")

        # u and fill should be renumberd in the same way, but types
        # must remain different, to let explicit u=0
        # self.values = map(lambda t: (f(t[0], t[1]), t[1]), self.values)
        newvals = []
        for t in self.values:
            if t[1] == "fill":
                t1 = "u"
            else:
                t1 = t[1]
            # newvals.append((f(t[0], t1), t[1]))
            if t1 in f:
                newval = f[t1](t[0])
            else:
                newval = t[0]
            newvals.append((newval, t[1]))
        self.values = newvals
        self.print_debug("after apply_map", "vi")
        return


# def _parse_geom(geom):
#     """
#     Parse the geometry part of a cell card.
#     """
#     raise NotImplementedError()
#     t = geom.split()
#     vals = []
#     fmts = []
#
#     # cell name
#     js = t.pop(0)
#     geom = geom.replace(js, tp, 1)
#     vals.append((int(js), 'cel'))
#     fmts.append(fmt_d(js))
#
#     if 'like' in geom.lower():
#         # this is like-but syntax
#         pass
#     else:
#         # get material and density.
#         # Density, if specified in cells card, should be allready hidden
#         ms = t.pop(0)
#         if int(ms) == 0:
#             inpt = inpt.replace(ms, tp+tp , 1)
#         else:
#             inpt = inpt.replace(ms, tp, 1)
#             inpt = inpt.replace('~', '~'+tp, 1)
#         vals.append((int(ms), 'mat'))
#         fmts.append(fmt_d(ms))
#
#         # placeholder for geometry prefix
#         vals.append(('', '#gpr'))
#         fmts.append('{}')


def _split_cell(input_, self):
    """
    Replace integers in the meaningful parts of a cell card with format
    specifiers, and return a list of replaced values together with their types.

    """

    # meaningful parts together. Originally, they cannot have \n chars, since
    # all of them should land to the card template, therefore, after all
    # entries are replaced with format specifiers, it can be split back to a
    # list easily at \n positions.
    inpt = "\n".join(input_)

    vals = []  # list of values
    fmts = []  # value format. It has digits, thus inserted into inpt later.
    tp = "_"  # temporary placeholder for format specifiers

    # Parse part before parameters. This is different for usual and like-but
    # syntax.  As result, all entries are replaced in inpt and i, the index
    # where parameter's part starts in inpt, is computed.

    if "like " in inpt.lower():

        # Get cell name
        t = inpt.split()
        js = t.pop(0)
        inpt = inpt.replace(js, tp, 1)
        vals.append((int(js), "cel"))
        fmts.append(fmt_d(js))

        # Get reference cell name:
        t.pop(0)  # like
        js = t.pop(0)
        inpt = inpt.replace(js, tp, 1)
        vals.append((int(js), "cel"))
        fmts.append(fmt_d(js))

        # compute i -- where first param token starts
        t.pop(0)  # but
        p0 = t.pop(0)
        i = inpt.index(p0)
        parm = [p0] + t

    else:
        # cell card has usual format.

        t = inpt.split()

        # Get cell name
        js = t.pop(0)
        inpt = inpt.replace(js, tp, 1)
        vals.append((int(js), "cel"))
        fmts.append(fmt_d(js))

        # get material and density.
        # Density, if specified in cells card, should be allready hidden
        ms = t.pop(0)
        if int(ms) == 0:
            inpt = inpt.replace(ms, tp + tp, 1)
        else:
            inpt = inpt.replace(ms, tp, 1)
            inpt = inpt.replace("~", "~" + tp, 1)
        vals.append((int(ms), "mat"))
        fmts.append(fmt_d(ms))

        # placeholder for geometry prefix
        vals.append(("", "#gpr"))
        fmts.append("{}")

        # Get geometry and parameters blocks. I assume that geom and param
        # blocks are separated by at least one space, so there will be an
        # element in t starting with alpha char -- This will be the first token
        # from the param block.
        geom = []
        parm = []
        while t:
            e = t.pop(0)
            if e[0].isalpha() or e[0] == "*":
                parm = [e] + t
                break
            else:
                geom.append(e)

        # print '_split_cell geom', geom, parm
        # replace integer entries in geom block:
        for s in re_int.findall(" ".join(geom)):
            # print 's from re_int', repr(s)
            # s is a surface or a cell (later only if prefixed by #)
            t = "cel" if s[0] == "#" else "sur"
            s = s if s[0].isdigit() else s[1:]
            f = fmt_d(s)
            inpt = inpt.replace(s, tp, 1)
            # print 't', repr(t)
            # print 's', repr(s)
            # print 'f', repr(f)
            # print repr(inpt)
            vals.append((int(s), t))
            fmts.append(f)

        # geometry suffix
        vals.append(("", "#gsu"))
        fmts.append("{}")
        # insert placeholder for geometry suffix
        if parm:
            inpt = inpt.replace(parm[0], "_" + parm[0], 1)

        # At this point all geom entries are replaced in inpt. The rest should
        # work only with the parm part of inpt. To ensure this, inpt is splitted
        # into inpt_geom and inpt_parm:
        if parm:
            i = inpt.index(parm[0])
        else:
            i = len(inpt)

    inpt_geom = inpt[:i]
    inpt_parm = inpt[i:]

    # print 'From parsing'
    # print repr(inpt)
    # print i
    # print repr(inpt_geom)
    # print repr(inpt_parm)
    # print parm

    # replace values in parameters block. Values are prefixed with = or space(s)
    # Note that tmp and imp values must be hidden
    t = " ".join(parm).replace("=", " ").split()  # get rid of =.
    while t:
        s = t.pop(0)
        # print '_split_cell s: ', repr(s)
        if s.lower() == "u":
            vs = t.pop(0)
            vv = int(vs)
            vf = fmt_d(vs)
            vt = "u"
            inpt_parm = inpt_parm.replace(vs, tp, 1)
            vals.append((vv, vt))
            fmts.append(vf)
        elif "fill" in s.lower():
            # print '_split_cell: has fill!'
            # assume that only one integer follows the fill keyword, optionally
            # with transformation in parentheses.

            # set transformation dimension deg if *FILL or cos if FILL
            if "*" in s:
                vv = "*"
            else:
                vv = ""
            vf = fmt_s(vv)
            vt = "#tunit"
            vals.append((vv, vt))
            fmts.append(vf)

            vs = t.pop(0)
            # if transformation in parentheses follows the universe number
            # immediately, split this manually:
            if "(" in vs:
                i = vs.index("(")
                ttt = vs[i:]
                vs = vs[:i]
                # vs, ttt = vs.split('(')
                t.insert(0, ttt)
            vv = int(vs)
            vf = fmt_d(vs)
            vt = "fill"
            inpt_parm = inpt_parm.replace(vs, tp, 1)
            vals.append((vv, vt))
            fmts.append(vf)
            # fill value can be followed by transformation in parentheses
            # Fill value can be optionally followed by transformation number of
            # transformation parameters in parentheses
            if t and "(" in t[0]:
                vsl = []  # lists of strings, values, formats and types
                vvl = []
                vfl = []
                vtl = []

                # add opening parenthesis
                vsl.append("(")
                vvl.append("(")
                vfl.append(fmt_s("("))
                vtl.append("#(")  # #-types are internal, don't output in --mode info.
                t[0] = t[0].replace("(", "", 1)

                # add entries in parentheses and the closing parenthis
                while vsl[-1] != ")":
                    vs = t.pop(0)
                    if ")" in vs:
                        vs = vs.replace(")", "", 1)
                        if vs:
                            vsl.append(vs)
                            vvl.append(vs)
                            vfl.append(fmt_s(vs))
                            vtl.append("#tparam")
                        vsl.append(")")
                        vvl.append(")")
                        vfl.append(fmt_s(")"))
                        vtl.append("#)")
                    elif vs:
                        vsl.append(vs)
                        vvl.append(vs)
                        vfl.append(fmt_s(vs))
                        vtl.append("#tparam")

                # check if only one parameter in parenthethes -- it is tr
                # number, not tr parameter
                if len(vsl) == 3:
                    vvl[1] = int(vvl[1])
                    vfl[1] = fmt_d(vsl[1])
                    vtl[1] = "tr"

                # add all strings, values, formats and types:
                for vs, vv, vf, vt in zip(vsl, vvl, vfl, vtl):
                    inpt_parm = inpt_parm.replace(vs, tp, 1)  # TODO: here only parm part of inpt should be modified.
                    vals.append((vv, vt))
                    fmts.append(vf)

            # warn if there is possibility for an array following the fill
            # keyword:
            # TODO fill value can be an array
            if "fill" == s.lower() and "lat" in "".join(parm).lower():
                print("WARNING: fill keyword followed by an array", end=" ")
                print("cannot be parsed")

    inpt = inpt_geom + inpt_parm

    # replace '_' with fmts:
    for f in fmts:
        inpt = inpt.replace(tp, f, 1)

    return inpt.split("\n"), vals


def _split_surface(input_):
    """
    Similar to _split_cell(), but for surface cards.
    """
    inpt = "\n".join(input_)
    t = inpt.split()

    vals = []  # like in split_cell()
    fmts = []
    tp = "_"

    # get surface name:
    js = t.pop(0)
    if not js[0].isdigit():
        js = js[1:]
    inpt = inpt.replace(js, tp, 1)
    vals.append((int(js), "sur"))
    fmts.append(fmt_d(js))

    # get TR or periodic surface:
    ns = t.pop(0)
    if ns[0].isdigit():
        # TR is given
        inpt = inpt.replace(ns, tp, 1)
        vals.append((int(ns), "tr"))
        fmts.append(fmt_d(ns))
        st = t.pop(0)
    elif ns[0] == "-":
        # periodic surface
        ns = ns[1:]
        inpt = inpt.replace(ns, tp, 1)
        vals.append((int(ns), "sur"))
        fmts.append(fmt_d(ns))
        st = t.pop(0)
    elif ns[0].isalpha():
        # ns is the surface type
        st = ns
    else:
        raise ValueError(input_, inpt, ns)

    # define coefficients
    scoef = list(map(float, t))

    for f in fmts:
        inpt = inpt.replace(tp, f, 1)

    return inpt.split("\n"), vals, st, scoef


def _get_int(s):
    r = ""
    for c in s:
        if r and c.isalpha():
            break
        elif c.isdigit():
            r += c
    return r


def _parse_tr(input_):
    """
    input_ should be already passed through _split_data()
    """
    inpt = "\n".join(input_)
    inp1, inp2 = inpt.split(None, 1)
    if "*" in inp1:
        unit = "*"
    else:
        unit = ""

    svals = inp2.split()
    for s in svals:
        inp2 = inp2.replace(s, "{}", 1)

    fvals = [(float(s), "float") for s in svals]
    return unit, (inp1 + " " + inp2).split("\n"), fvals


def _split_data(input_):
    inpt = "\n".join(input_)
    t = inpt.split()

    vals = []
    fmts = []
    tp = "_"

    if "tr" in t[0][:3].lower():
        # TRn card
        dtype = "TRn"
        ns = _get_int(t[0])
        inpt = inpt.replace(ns, tp, 1)
        vals.append((int(ns), "tr"))
        fmts.append(fmt_d(ns))
    elif t[0][0].lower() == "m" and "mode" not in t[0].lower() and "mesh" not in t[0].lower() and "mphys" not in t[0].lower():
        # This is the Mn, MTn or MPNn card
        ms = _get_int(t[0])
        inpt = inpt.replace(ms, tp, 1)
        vals.append((int(ms), "mat"))
        fmts.append(fmt_d(ms))
        # additional tests to define data card type:
        if t[0][1].isdigit():
            dtype = "Mn"
        elif t[0][1].lower() == "t":
            dtype = "MTn"
        elif t[0][1].lower() == "p":
            dtype = "MPNn"
        elif t[0][1].lower() == "x":
            dtype = "MXn"
        else:
            dtype = None
    elif t[0][0].lower() == "f" and t[0][1].isdigit():
        # FN card
        dtype = "Fn"
        ns = _get_int(t[0])  # tally number
        inpt = inpt.replace(ns, tp, 1)
        vals.append((int(ns), "tal"))
        fmts.append(fmt_d(ns))

        # define type of integers by tally type:
        nv = int(ns[-1])
        if nv in [1, 2]:
            typ = "sur"
        elif nv in [4, 6, 7, 8]:
            typ = "cel"
        else:
            typ = ""

        if typ:
            # Lattice indices, surrounded by square brakets must allready be
            # hidden

            # Special treatment, if tally has 'u=' syntax.
            hasu = "u" in inpt.lower() and "=" in inpt.lower()
            # find all integers -- they are cells or surfaces
            for s in re_int.findall(inpt):
                ss = s[1:]
                tpe = typ
                if hasu:
                    # ss can be universe. To distinguish this, one needs to look
                    # back in previous cheracters in c.
                    i1 = inpt.rfind(tp)
                    i2 = inpt.find(ss)
                    part = inpt[i1:i2]
                    while " " in part:
                        part = part.replace(" ", "")
                    if part[-2:].lower() == "u=":
                        tpe = "u"
                inpt = inpt.replace(ss, tp, 1)
                vals.append((int(ss), tpe))
                fmts.append(fmt_d(ss))
    elif "fmesh" == t[0][:5].lower() and t[0][5].isdigit():
        # fmesh card
        dtype = "fmesh"
        ns = _get_int(t[0])  # tally number
        inpt = inpt.replace(ns, tp, 1)
        vals.append((int(ns), "tal"))
        fmts.append(fmt_d(ns))
    else:
        dtype = None

    for f in fmts:
        inpt = inpt.replace(tp, f, 1)

    return inpt.split("\n"), vals, dtype


def is_commented(l):
    """
    Return True if l is a commented line.
    """
    res = False

    # remove newline chars at the end of l:
    l = l.splitlines()[0]
    ls = l[0:6].lstrip().lower()
    if "c " == ls[0:2]:
        res = True
        # print 'is_com "c "',
    elif "c" == l.lower():
        res = True
        # print 'is_com "c"',
    # print 'is_com', res
    return res


def is_fc_card(l):
    """
    Return true, if line l is tally comment cards, fcN
    """
    return l.lstrip().lower()[:2] == "fc"


def is_blankline(l):
    """
    Return True, if l is the delimiter blank line.
    """
    return l.strip() == ""


def get_cards(inp, debug=None, preservetabs=False):
    """
    Check first existence of a dump file

    If dump exists and it is newwer than the input file, read the dump file
    """
    for c in get_cards_from_input(inp, debug=debug, preservetabs=preservetabs):
        yield c


def index_(line, chars="$&"):
    """
    Find the first index of one of the chars in line.
    """
    r = re.compile(f"[{chars}]")
    m = r.search(line)
    if m:
        i = m.end() - 1
    else:
        i = len(line) - 1
    return i


def get_cards_from_input(inp, debug=None, preservetabs=False):
    """
    Iterable, return instances of the Card() class representing
    cards in the input file.

    inp -- is the filename.
    """

    def _yield(card, ct, ln):
        return Card(card, ct, ln, debug)

    def replace_tab(l, cln, preserve=False, ts=8):
        """
        Replace tabs as in MCNP5 (Vol II, Chapter 1 - Primer, I. MCNP INPUT FOR
        SAMPLE PROBLEM, A. INP File, p. 1-3)
        """
        if preserve:
            return l[:]
        else:
            while "\t" in l:
                i = l.index("\t")
                ii = (i // ts + 1) * ts - i
                print(f"c Line {cln + 1}: tab replaced with {ii} spaces")
                l = l[:i] + " " * ii + l[i + 1 :]
            return l[:]

    cln = 0  # current line number. Used only for debug
    with open(inp, "r") as f:
        # define the first block:
        # -----------------------

        # Next block ID
        ncid = 0  # 0 is not used in card ID dictionary CID.

        # Parse the 1-st line. It can be message, cell or data block.
        l = replace_tab(next(f), cln, preserve=preservetabs)
        cln += 1
        # kw = l.lower().split()[0]
        kw = l.lstrip()
        if "message:" == kw[:8].lower():
            # read message block right here
            res = []
            while not is_blankline(l):
                res.append(l)
                l = replace_tab(next(f), cln, preserve=preservetabs)
                cln += 1
            yield _yield(res, CID.message, cln - 1)  # message card
            yield _yield(l, CID.blankline, cln)  # blank line
            l = replace_tab(next(f), cln, preserve=preservetabs)
            cln += 1
            ncid = CID.title
        elif "continue" == kw[:8].lower():
            # input file for continue job. Contains only data block.
            ncid = CID.data
        else:
            ncid = CID.title
        if ncid == CID.title:
            # l contains the title card
            yield _yield([l], ncid, cln)
            ncid += 1

        # read all other lines
        # --------------------

        # Line can be a continuation line in the following cases:
        #   * all lines in the message block, i.e. before the blank line
        #     delimiter
        #   * if line starts with 5 or more spaces,
        #   * if previous line ends with & sign.
        # Thus, the role of the current line (continuation or not) can be
        # defined by the line itself (5 spaces), or by previous lines (message
        # block or & sign). This can lead to inconsistency, when previous line
        # is delimited by &, but is followed by the blank line delimiter.  in
        # this case (rather theretical), blank line delimiter (as appeared more
        # lately) delimites the card from the previous line.
        cf = False  # continuation line flag. True only when prev. line has &.

        # Comment lines (CL) can be between cards or inside them. CL between two
        # cards are yielded as block of comments (although usually, CL are used
        # to describe cards that follow them).  CL inside a card will belong to
        # the card.

        card = []  # card is a list of lines.
        cmnt = []  # list of comment lines.
        for l in f:
            l = replace_tab(l, cln, preserve=preservetabs)
            cln += 1
            if is_blankline(l):
                # blank line delimiter. Stops card even if previous line
                # contains &
                if card:
                    # card can be empty, for example, when several empty lines
                    # are at the end of file
                    yield _yield(card, ncid, cln - len(card) - len(cmnt))
                if cmnt:
                    yield _yield(cmnt, CID.comment, cln - len(cmnt))
                    cmnt = []
                yield _yield(l, CID.blankline, cln)
                ncid += 1
                card = []
                if ncid == 6:
                    break
            elif l[0:5] == "     " or cf:
                # l is continuation line.
                if cmnt:
                    card += cmnt  # prev. comment lines belong to this card.
                    cmnt = []
                card.append(l)
                cf = l[: index_(l)].find("&", 0, 81) > -1
            elif is_commented(l):
                # l is a line comment. Where it belongs (to the current card or
                # to the next one), depends on the next line, therefore, just
                # store temorarily.
                cmnt.append(l)
            else:
                # l is the 1-st line of a card. Return previous card and
                # comments
                if card:
                    yield _yield(card, ncid, cln - len(card) - len(cmnt))
                if cmnt:
                    yield _yield(cmnt, CID.comment, cln - len(cmnt))
                    cmnt = []
                card = [l]
                # if tally comment card, i.e. started with fc, the & character
                # does not mean continuation.
                cf = not is_fc_card(l) and l[: index_(l)].find("&", 0, 81) > -1
        if card:
            yield _yield(card, ncid, cln - len(card) - len(cmnt))
        if cmnt:
            yield _yield(cmnt, CID.comment, cln - len(cmnt))


def get_blocks(cards):
    """
    Return a dict of cards in blocks.
    """

    d = {}
    cbt = None  # current block type
    cbc = []  # current block cards
    for c in cards:
        if c.ctype == CID.blankline:
            d[cbt] = cbc
            cbt = None
            cbc = []
        elif c.ctype == CID.title:
            d[c.ctype] = [c]
        else:
            cbc.append(c)
            if cbt is None and c.ctype > 0:
                cbt = c.ctype
    if cbc:
        d[cbt] = cbc
    return d


def are_close_vals(x, y, re=1e-6, ra=0.0):
    """
    Return True if x and y are closer then re or ra.
    """
    if abs(x - y) <= ra:
        r = True
    elif x != 0:
        r = abs((x - y) / x) <= re
    else:
        # y is not equal to x and x is 0 -> y is not 0.
        r = abs((x - y) / y) <= re
    return r


def are_close_lists(x, y, re=1e-6, pci=[]):
    """
    Return True if x and y are close but not equal.
    """
    if len(x) != len(y):
        res = False
        msg = "Different length"

    if x == y:
        return True

    # pci -- list of indices that define elements of x and y to be checked for
    # proportionality only.
    if len(pci) == 0:
        # empty list means all x and y elements compare without arbitrary
        # normalization.
        xe = x[:]
        ye = y[:]
        xp = []
        yp = []
    else:
        if len(pci) % 2 == 1:
            # augment with len(x) +1
            pci = tuple(pci) + (len(x) + 1,)
        xe = []
        ye = []
        xp = []
        yp = []
        i = 0
        for i1, i2 in zip(pci[0::2], pci[1::2]):
            xe += x[i:i1]
            ye += y[i:i1]
            xp += x[i1:i2]
            yp += y[i1:i2]
            i = i2

    # normalize yp
    xpn = sum([e**2 for e in xp])
    ypn = sum([e**2 for e in yp])
    if xpn > 0 and ypn > 0:
        yp = [e * xpn / ypn for e in yp]

    msg = []
    res = []
    for xl, yl in zip([xe, xp], [ye, yp]):
        # compare xl and yl without normalization
        if xl == yl:
            res.append(True)
            msg.append("exact match")
        else:
            n = 0
            for xx, yy in zip(xl, yl):
                r = are_close_vals(xx, yy, re)
                if not r:
                    m = f"diff at {n}"
                    break
            else:
                m = "all elements are close or equal"
                r = True
            res.append(r)
            msg.append(m)

        if not res[-1]:
            result = False
            break

    else:
        result = True
    return result


if __name__ == "__main__":
    pass
