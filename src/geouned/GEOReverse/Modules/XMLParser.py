from .remh import Cline


class CellCard:

    def __init__(self, data):
        self.type = "cell"
        self.TR = None
        self.processData(data)

    def processData(self, data):
        self.name = int(data["id"])
        self.level = None

        if "material" in data.keys():
            self.MAT = 0 if data["material"] == "void" else int(data["material"])
        else:
            self.MAT = None

        if "universe" in data.keys():
            self.U = 0 if int(data["universe"]) == 1 else int(data["universe"])
        else:
            self.U = 0

        if "fill" in data.keys():
            self.FILL = int(data["fill"])
        else:
            self.FILL = None

        self.geom = Cline(data["region"].replace("|", ":"))


class SurfCard:
    def __init__(self, data):

        self.type = "surface"
        self.processData(data)

    def processData(self, data):
        self.name = int(data["id"])
        self.stype = data["type"]
        self.scoefs = tuple(float(x) for x in data["coeffs"].split())


def get_cards(root):
    for c in root:
        yield process_card(c)


def process_card(card):
    ctype = card.tag
    if ctype == "cell":
        return CellCard(card.attrib)

    elif ctype == "surface":
        return SurfCard(card.attrib)
