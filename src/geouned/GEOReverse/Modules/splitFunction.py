import math

import BOPTools.SplitAPI
import FreeCAD
import Part


class SplitBase:
    def __init__(self, base, knownSurf={}):
        self.base = base
        self.knownSurf = knownSurf


def joinBase(baseList):
    shape = []
    surf = {}
    removedKeys = []
    for b in baseList:
        if b.base is not None:
            shape.append(b.base)
        for k, v in b.knownSurf.items():
            if k in removedKeys:
                continue
            if k not in surf.keys():
                surf[k] = v
            else:
                if surf[k] == v:
                    continue
                else:
                    surf[k] = None
                    removedKeys.append(k)

    newbase = FuseSolid(shape)
    return SplitBase(newbase, surf)


# TODO rename this function as there are two with the name name
def SplitSolid(base, surfacesCut, cellObj, solidTool=False, tolerance=0.01):  # 1e-2
    # split Base (shape Object or list/tuple of shapes)
    # with selected surfaces (list of surfaces objects) cutting the base(s) (surfacesCut)
    # cellObj is the CAD object of the working cell to reconstruction.
    # the function return a list of solids enclosed fully in the cell (fullPart)
    # and a list of solids not fully enclosed in the cell (cutPart). These lasts
    # will require more splitting with the others surfaces defining the cell.

    fullPart = []
    cutPart = []

    # part if several base in input
    if type(base) is list or type(base) is tuple:
        for b in base:
            fullList, cutList = SplitSolid(b, surfacesCut, cellObj, tolerance=tolerance)
            fullPart.extend(fullList)
            cutPart.extend(cutList)
        return fullPart, cutPart

    # part if base is shape object

    if solidTool:
        Tools = (cellObj.shape,)
    else:
        Tools = tuple(s.shape for s in surfacesCut)
    # for s in surfacesCut:
    #    print(s.type,s.params,s.id)
    #    s.shape.exportStep('tool{}.stp'.format(s.id))
    # base.base.exportStep('base.stp')
    Solids = BOPTools.SplitAPI.slice(base.base, Tools, "Split", tolerance=tolerance).Solids
    if not Solids:
        Solids = [base.base]
    partPositions, partSolids = space_decomposition(Solids, surfacesCut)

    for pos, sol in zip(partPositions, partSolids):
        # fullPos = updateSurfacesValues(pos,cellObj.surfaces,base.knownSurf)
        # inSolid = cellObj.definition.evaluate(fullPos)

        if not solidTool:
            pos.update(base.knownSurf)
        inSolid = cellObj.definition.evaluate(pos)

        # if solidTool :
        #  ii += 1
        #  print(solidTool)
        #  print(cellObj.definition)
        #  print(pos)
        #  print('eval',inSolid)
        #  name = str(cellObj.definition)
        #  sol.exportStep('solid_{}{}.stp'.format(name,ii))

        if inSolid:
            fullPart.append(SplitBase(sol, pos))
        elif inSolid is None:
            cutPart.append(SplitBase(sol, pos))
    return fullPart, cutPart


def updateSurfacesValues(position, surfaces, knownSurf):
    position.update(knownSurf)
    sname = set(surfaces.keys())
    pname = set(position.keys())

    fullpos = position.copy()
    for name in sname.difference(pname):
        fullpos[name] = None
    return fullpos


# Get the position of subregion with respect
# all cutting surfaces
def space_decomposition(solids, surfaces):

    component = []
    good_solids = []
    for c in solids:
        if c.Volume < 1e-3:
            if abs(c.Volume) < 1e-3:
                continue
            else:
                c.reverse()
                print("Negative solid Volume", c.Volume)
        Svalues = {}
        point = point_inside(c)
        if point == None:
            continue  # point not found in solid (solid is surface or very thin can be source of lost particules in MCNP)
        for surf in surfaces:
            Svalues[surf.id] = surface_side(point, surf)

        component.append(Svalues)
        good_solids.append(c)
    return component, good_solids


# find one point inside a solid (region)
def point_inside(solid):

    cut_line = 32
    cut_box = 4

    # no poner boundbox, el punto puente caer en una superficie para geometria triangular
    point = solid.CenterOfMass
    if solid.isInside(point, 0.0, False):
        return point

    v1 = solid.Vertexes[0].Point
    for vi in range(len(solid.Vertexes) - 1, 0, -1):
        v2 = solid.Vertexes[vi].Point
        dv = (v2 - v1) * 0.5

        n = 1
        while True:
            for i in range(n):
                point = v1 + dv * (1 + 0.5 * i)
                if solid.isInside(point, 0.0, False):
                    return point
            n = n * 2
            dv = dv * 0.5
            if n > cut_line:
                break

    #      Box_Volume = BBox.XLength*BBox.YLength*BBox.ZLength
    #      if (solid.Volume < Box_Volume/ math.pow(16,nmax_cut)) :
    #           print('very small Solid Volume (solid volume, box volume): {},{}'.format(solid.Volume,Box_Volume))
    #           return None
    BBox = solid.optimalBoundingBox(False)
    box = [BBox.XMin, BBox.XMax, BBox.YMin, BBox.YMax, BBox.ZMin, BBox.ZMax]

    boxes, centers = divide_box(box)
    n = 0

    while True:
        for p in centers:
            pp = FreeCAD.Vector(p[0], p[1], p[2])
            if solid.isInside(pp, 0.0, False):
                return pp

        subbox = []
        centers = []
        for b in boxes:
            btab, ctab = divide_box(b)
            subbox.extend(btab)
            centers.extend(ctab)
        boxes = subbox
        n = n + 1

        if n == cut_box:
            print(f"Solid not found in bounding Box (Volume : {solid.Volume})")
            print("Valid Solid : ", solid.isValid())
            return None


# divide a box into 8 smaller boxes
def divide_box(Box):
    xmid = (Box[1] + Box[0]) * 0.5
    ymid = (Box[3] + Box[2]) * 0.5
    zmid = (Box[5] + Box[4]) * 0.5

    b1 = (Box[0], xmid, Box[2], ymid, Box[4], zmid)
    p1 = (0.5 * (Box[0] + xmid), 0.5 * (Box[2] + ymid), 0.5 * (Box[4] + zmid))

    b2 = (xmid, Box[1], Box[2], ymid, Box[4], zmid)
    p2 = (0.5 * (xmid + Box[1]), 0.5 * (Box[2] + ymid), 0.5 * (Box[4] + zmid))

    b3 = (Box[0], xmid, ymid, Box[3], Box[4], zmid)
    p3 = (0.5 * (Box[0] + xmid), 0.5 * (ymid + Box[3]), 0.5 * (Box[4] + zmid))

    b4 = (xmid, Box[1], ymid, Box[3], Box[4], zmid)
    p4 = (0.5 * (xmid + Box[1]), 0.5 * (ymid + Box[3]), 0.5 * (Box[4] + zmid))

    b5 = (Box[0], xmid, Box[2], ymid, zmid, Box[5])
    p5 = (0.5 * (Box[0] + xmid), 0.5 * (Box[2] + ymid), 0.5 * (zmid + Box[5]))

    b6 = (xmid, Box[1], Box[2], ymid, zmid, Box[5])
    p6 = (0.5 * (xmid + Box[1]), 0.5 * (Box[2] + ymid), 0.5 * (zmid + Box[5]))

    b7 = (Box[0], xmid, ymid, Box[3], zmid, Box[5])
    p7 = (0.5 * (Box[0] + xmid), 0.5 * (ymid + Box[3]), 0.5 * (zmid + Box[5]))

    b8 = (xmid, Box[1], ymid, Box[3], zmid, Box[5])
    p8 = (0.5 * (xmid + Box[1]), 0.5 * (ymid + Box[3]), 0.5 * (zmid + Box[5]))

    return [b1, b2, b3, b4, b5, b6, b7, b8], [p1, p2, p3, p4, p5, p6, p7, p8]


# check the position of the point with respect
# a surface
def surface_side(p, surf):
    if surf.type == "sphere":
        org, R = surf.params
        D = p - org
        inout = D.Length - R

    elif surf.type == "plane":
        normal, d = surf.params
        inout = p.dot(normal) - d

    elif surf.type == "cylinder":
        P, v, R = surf.params

        D = p - P
        if not surf.truncated:
            inout = D.cross(v).Length - R
        else:
            inCyl = D.cross(v).Length / v.Length - R  # <0 in cylinder
            inPln = btwPPlanes(p, P, v)  # <0  between planes

            if (inCyl < 0) and (inPln < 0):
                inout = -1  # inside the can
            else:
                inout = 1  # outside the can

    elif surf.type == "cone":
        if not surf.truncated:
            P, v, t, dblsht = surf.params
            X = p - P
            X.normalize()
            dprod = X.dot(v)
            dprod = max(-1, min(1, dprod))
            a = math.acos(dprod) if not dblsht else math.acos(abs(dprod))
            inout = a - math.atan(t)
        else:
            P, v, R1, R2 = surf.params
            apex = P + R1 / (R1 - R2) * v

            X = p - apex
            X.normalize()
            dprod = X.dot(-v) / v.Length  # -v because reverse axis. in MCNP TRC r1 > r2
            dprod = max(-1, min(1, dprod))
            a = math.acos(dprod)

            t = (R1 - R2) / v.Length
            inCone = a - math.atan(t)
            inPln = btwPPlanes(p, P, v)  # <0  between planes

            if (inCone < 0) and (inPln < 0):
                inout = -1  # inside the can
            else:
                inout = 1  # outside the can

    elif surf.type == "cone_elliptic":
        apex, axis, Ra, radii, rAxes, dblsht = surf.params

        r = p - apex
        X = r.dot(rAxes[1])
        Y = r.dot(rAxes[0])
        Z = r.dot(axis)
        if dblsht:
            Z = abs(Z)
        inout = (X / radii[1]) ** 2 + (Y / radii[0]) ** 2 - Z / Ra

    elif surf.type == "hyperboloid":
        center, axis, radii, rAxes, onesht = surf.params

        r = p - center
        rX = r.dot(rAxes[1])
        v = r - (rX * rAxes[1] + center)
        d = v.Length

        one = 1 if onesht else -1
        radical = (rX / radii[1]) ** 2 + one

        if radical > 0:
            Y = radii[0] * math.sqrt(radical)
            inout = d - Y
        else:
            inout = 1

    elif surf.type == "ellipsoid":
        center, axis, radii, rAxes = surf.params

        r = p - center
        rX = r.dot(axis)
        rY = r - (rX * axis + center)

        if axis.add(-rAxes[0]).Length < 1e-5:
            radX, radY = radii
        else:
            radY, radY = radii

        radical = 1 - (rX / radX) ** 2
        if radical > 0:
            Y = radY * math.sqrt(radical)
            inout = rY - Y
        else:
            inout = 1

    elif surf.type == "cylinder_elliptic":
        center, axis, radii, rAxes = surf.params

        r = p - center
        X = r.dot(rAxes[1])
        Y = r.dot(rAxes[0])
        inout = (X / radii[1]) ** 2 + (Y / radii[0]) ** 2 - 1

        if surf.truncated and inout < 0:
            inout = btwPPlanes(p, center, axis)  # <0  between planes

    elif surf.type == "cylinder_hyperbolic":
        center, axis, radii, rAxes = surf.params

        r = p - center
        X = r.dot(rAxes[1])
        Y = r.dot(rAxes[0])
        inout = (X / radii[1]) ** 2 - (Y / radii[0]) ** 2 - 1

    elif surf.type == "paraboloid":
        center, axis, focal = surf.params

        r = p - center
        X = r.dot(axis)
        if X < 0:
            inout = 1
        else:
            v = r - X * axis
            d = v.Length
            Y = math.sqrt(4 * focal * X)
            inout = d - Y

    elif surf.type == "torus":
        P, v, Ra, Rb, Rc = surf.params

        d = p - P
        z = d.dot(v)
        rz = d - z * v
        inout = (z / Rb) ** 2 + ((rz.Length - Ra) / Rc) ** 2 - 1

    elif surf.type == "box":
        P, v1, v2, v3 = surf.params
        for v in (v1, v2, v3):
            inout = btwPPlanes(p, P, v)  # <0  between planes
            if inout > 0:
                break

    else:
        print(f"surface type {surf[0]} not considered")
        return

    return inout > 0


def btwPPlanes(p, p0, v):

    p1 = p0 + v
    inP0 = v.dot(p - p0)  # >0 plane(base plane) side inside the cylinder
    inP1 = v.dot(p - p1)  # >0 plane(base plane) side inside the cylinder

    if (inP0 > 0) and (inP1 < 0):
        return -1
    else:
        return 1


# ************************************************


def FuseSolid(parts):
    if (len(parts)) <= 1:
        if parts:
            solid = parts[0]
        else:
            return None
    else:
        try:
            fused = parts[0].fuse(parts[1:])
        except:
            fused = None

        if fused is not None:
            try:
                refinedfused = fused.removeSplitter()
            except:
                refinedfused = fused

            if refinedfused.isValid():
                solid = refinedfused
            else:
                if fused.isValid():
                    solid = fused
                else:
                    solid = Part.makeCompound(parts)
        else:
            solid = Part.makeCompound(parts)

    if solid.Volume < 0:
        solid.reverse()
    return solid
