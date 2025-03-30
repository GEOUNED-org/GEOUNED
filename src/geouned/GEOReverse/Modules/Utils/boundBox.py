import FreeCAD
import Part
import math
import numpy
from .booleanFunction import BoolSequence
from ..data_class import BoxSettings

twoPi = math.pi * 2


class myBox:
    def __init__(self, boundBox=None, orientation=None):

        if boundBox is not None:
            if boundBox.XLength <= 1e-12:
                self.Box = None
            elif boundBox.YLength <= 1e-12:
                self.Box = None
            elif boundBox.ZLength <= 1e-12:
                self.Box = None
            else:
                self.Box = boundBox
        else:
            self.Box = None
        self.Orientation = orientation

    def add(self, box):
        if self.Orientation is None:
            self.Box = box.Box
            self.Orientation = box.Orientation
        elif self.Box is None:
            if self.Orientation == "Forward":
                self.Box = box.Box
                self.Orientation = box.Orientation
        elif box.Box is None:
            if box.Orientation == "Reversed":
                self.Box = None
                self.Orientation = "Reversed"
        elif self.Orientation == box.Orientation:
            self.Box.add(box.Box)
        else:
            # -A OR B == -(A AND -B)
            if self.Orientation == "Forward":
                Rbox, Fbox = self, box
            else:
                Rbox, Fbox = box, self
            self.Box = box_intersect(Fbox, Rbox)
            self.Orientation = "Reversed"

    def mult(self, box):
        if self.Orientation is None:
            self.Box = box.Box
            self.Orientation = box.Orientation
        elif self.Box is None:
            if self.Orientation == "Reversed":
                self.Box = box.Box
                self.Orientation = box.Orientation
        elif box.Box is None:
            if box.Orientation == "Forward":
                self.Box = None
                self.Orientation = "Forward"
        elif self.Orientation == box.Orientation:
            inter = self.Box.intersected(box.Box)
            if inter.isValid():
                self.Box = inter
            else:
                self.Box = None
        else:
            if self.Orientation == "Forward":
                Fbox, Rbox = self, box
            else:
                Fbox, Rbox = box, self
            self.Box = box_intersect(Fbox, Rbox)
            self.Orientation = "Forward"

    def sameBox(self, box):
        if self.Box is None or box.Box is None:
            if self.Box is None and box.Box is None:
                return self.Orientation == box.Orientation
            else:
                return False

        for i in range(6):
            p1 = self.Box.getPoint(i)
            p2 = box.Box.getPoint(i)
            if (p1 - p2).Length > 1e-6:
                return False
        return True


class solid_plane_box:
    def __init__(self, NTCell=None, outbox=None, orientation="Undefined"):
        if NTCell is None:
            self.planes = None
            self.definition = None
            self.surfaces = None
            self.surf_to_plane = None
            self.insolid_tolerance = BoxSettings().insolid_tolerance
            self.universe_radius = BoxSettings().universe_radius
            self.universe_center = FreeCAD.Vector(0, 0, 0)
            self.orientation = None
        else:
            self.insolid_tolerance = NTCell.settings.insolid_tolerance
            self.universe_radius = NTCell.settings.universe_radius
            self.surfaces = NTCell.surfaces
            plane_dict, surf_to_plane_dict = quadric_to_plane(NTCell.definition, NTCell.surfaces, orientation)
            self.planes = plane_dict
            self.surf_to_plane = surf_to_plane_dict
            self.definition = plane_definition(NTCell.definition.copy(), surf_to_plane_dict, orientation)
            self.orientation = self.get_box_orientation()

            if orientation != self.orientation:
                plane_dict, surf_to_plane_dict = quadric_to_plane(NTCell.definition, NTCell.surfaces, self.orientation)
                self.planes = plane_dict
                self.surf_to_plane = surf_to_plane_dict

            if orientation == "Undefined" and self.orientation != "Undefined":
                self.definition = plane_definition(NTCell.definition.copy(), surf_to_plane_dict, self.orientation)

            self.universe_center = FreeCAD.Vector(0, 0, 0)

        self.outBox = outbox
        if outbox:
            self.universe_radius = outbox.Box.DiagonalLength * 0.5
            self.universe_center = outbox.Box.Center
        else:
            r = self.universe_radius
            self.outBox = myBox(FreeCAD.BoundBox(-r, -r, -r, r, r, r), "Forward")

    def export_surf_planes(self, box):

        surf = set(self.surf_to_plane.keys())
        surf_planes = set()
        for s in surf:
            planes = []
            for p in self.surf_to_plane[s]:
                if p not in self.surf_to_plane.keys():
                    break
                normal, position = self.planes[p].Axis, self.planes[p].Position
                planes.append(makePlane(normal, position, box))
                surf_planes.add(p)
            else:
                compsurf = Part.Compound(planes)
                if compsurf is not None:
                    compsurf.exportStep(f"psurf_{s}.stp")

        all_planes = set(self.planes.keys())
        for p in all_planes - surf_planes:
            normal, position = self.planes[p].Axis, self.planes[p].Position
            compsurf = makePlane(normal, position, box)
            if compsurf is not None:
                compsurf.exportStep(f"psurf_{p}.stp")

    def isInside(self, point, boundary_inside=True):
        surf_value = dict()
        for p_index, p in self.planes.items():
            normal, pointPlane = p.Axis, p.Position
            pt = FreeCAD.Vector(point.x, point.y, point.z)
            r = pt - pointPlane
            dot = normal.dot(r)
            if abs(dot) < self.insolid_tolerance:
                surf_value[p_index] = None  # undefined value for point close to the surface
            else:
                surf_value[p_index] = dot > 0
        inside = self.definition.evaluate(surf_value)

        # if point close to the surface assume inside the solid independently if inside or outside
        # inside if point close to boundary
        forward_inside = boundary_inside if inside is None else inside
        return forward_inside if boundary_inside else not forward_inside

    def copy(self, newdefinition=None, rebuild=False):
        cpsol = solid_plane_box()
        cpsol.insolid_tolerance = self.insolid_tolerance
        cpsol.universe_radius = self.universe_radius
        cpsol.universe_center = self.universe_center
        cpsol.outBox = self.outBox
        cpsol.orientation = self.orientation

        if newdefinition is not None:
            if rebuild:
                plane_dict, surf_to_plane_dict = quadric_to_plane(newdefinition, self.surfaces, self.orientation)
                cpsol.planes = plane_dict
                cpsol.surf_to_plane = surf_to_plane_dict
                cpsol.definition = plane_definition(newdefinition.copy(), surf_to_plane_dict, self.orientation)
            else:
                cpsol.surf_to_plane = self.surf_to_plane
                cpsol.definition = newdefinition.copy()
                cpsol.planes = dict()
                for p in newdefinition.get_surfaces_numbers():
                    if p in self.planes.keys():
                        cpsol.planes[p] = self.planes[p]
        else:
            cpsol.surf_to_plane = self.surf_to_plane
            cpsol.definition = self.definition.copy()
            cpsol.planes = self.planes
        return cpsol

    def get_boundBox(self, enlarge=0):
        mBox = self.build_box_depth()
        bBox = mBox.Box
        if bBox is not None and enlarge > 0:
            dx = (bBox.XMax - bBox.XMin) * enlarge
            dy = (bBox.YMax - bBox.YMin) * enlarge
            dz = (bBox.ZMax - bBox.ZMin) * enlarge
            bBox = FreeCAD.BoundBox(
                bBox.XMin - dx, bBox.YMin - dy, bBox.ZMin - dz, bBox.XMax + dx, bBox.YMax + dy, bBox.ZMax + dz
            )
            mBox.Box = bBox
        return mBox

    def build_box_depth(self):

        if self.definition.level == 0:
            return self.get_component_boundBox()
        else:
            box_list = []
            if type(self.definition.elements) is bool:
                return myBox(None, "Reversed") if self.definition.elements else myBox(None, "Forward")
            for c in self.definition.elements:
                cbox = self.copy(c)
                box = cbox.build_box_depth()
                box_list.append(box)

            fullBox = myBox()
            if self.definition.operator == "AND":
                for box in box_list:
                    fullBox.mult(box)
                    if fullBox.Box is None and fullBox.Orientation == "Forward":
                        break
            else:
                for box in box_list:
                    fullBox.add(box)
                    if fullBox.Box is None and fullBox.Orientation == "Reversed":
                        break

            return fullBox

    def get_component_boundBox(self, cutBoundary=False):
        axis_list = ("x", "y", "z")

        orientation = self.get_box_orientation()
        if not cutBoundary:
            if orientation == "Undefined":
                cutBoundary = True
                orientation = "Forward"
        else:
            if orientation == "Undefined":
                orientation = "Forward"
        # if orientation == "Undefined":
        #    orientation = "Forward"

        planes_inter = tuple(self.planes[x] for x in self.definition.get_surfaces_numbers())
        point_list = plane_intersect(planes_inter, self.outBox.Box, cutBoundary)
        box_lim = []
        if len(point_list) < 6:
            if not cutBoundary:
                return self.get_component_boundBox(True)
            else:
                return myBox(None, orientation)
        else:
            for axis in axis_list:
                s_point = sort_point(point_list, axis)
                if s_point == []:
                    return None
                for point in s_point:
                    if self.isInside(point, orientation == "Forward"):
                        box_lim.append(pointaxis(point, axis))
                        break

                s_point = remove_points(s_point, pointaxis(point, axis), axis, True)
                if s_point == []:
                    return None
                for point in s_point:
                    if self.isInside(point, orientation == "Forward"):
                        box_lim.append(pointaxis(point, axis))
                        break
                point_list = remove_points(s_point, pointaxis(point, axis), axis, False)

            # if len(box_lim) < 6:
            #    return myBox(None, orientation)
            # else:
            #    box = FreeCAD.BoundBox(box_lim[0], box_lim[2], box_lim[4], box_lim[1], box_lim[3], box_lim[5])
            #    return myBox(box, orientation)

            if len(box_lim) < 6:
                if cutBoundary:
                    return myBox(None, orientation)
                else:
                    return self.get_component_boundBox(True)
            else:
                box = FreeCAD.BoundBox(box_lim[0], box_lim[2], box_lim[4], box_lim[1], box_lim[3], box_lim[5])
                if box.XLength < 1e-12 or box.YLength < 1e-12 or box.ZLength < 1e-12:
                    if cutBoundary:
                        return myBox(None, orientation)
                    else:
                        return self.get_component_boundBox(True)
                else:
                    return myBox(box, orientation)

    def get_box_orientation(self):
        ninside = 0
        universeBox = FreeCAD.BoundBox(
            -self.universe_radius,
            -self.universe_radius,
            -self.universe_radius,
            self.universe_radius,
            self.universe_radius,
            self.universe_radius,
        )
        for i in range(8):
            p = universeBox.getPoint(i)
            if self.isInside(p, True):
                ninside += 1
        if ninside == 8:
            return "Reversed"
        elif ninside == 0:
            return "Forward"
        else:
            return "Undefined"


def quadric_to_plane(cellDef, surfaces, orientation):

    surf_planes_dict = dict()
    planes = dict()

    surf_index = cellDef.signedSurfaces()
    next = list({abs(s) for s in surf_index})
    next.sort()
    next_index = next[-1] + 1
    apex = []

    if orientation == "Reversed":
        chg = True
    elif orientation == "Forward":
        chg = False
    else:
        chg = None

    for s_index in surf_index:
        s_index = abs(s_index)
        s = surfaces[s_index]
        if s.type == "plane":
            normal, d = s.params
            position = normal * d
            planes[s_index] = Part.Plane(position, normal)
        else:
            if chg is None:
                pos = None
            else:
                pos = (s_index > 0) == chg
            surf_planes = convert_to_planes(s, pos)
            if s.type == "cone":
                apex.append(s.params[0])
                dbl = s.params[3]
                p_index = []
                for p in surf_planes:
                    planes[next_index] = p
                    p_index.append(next_index)
                    next_index += 1

                if dbl:
                    surf_planes_dict[s_index] = ("dblcone", p_index)
                else:
                    surf_planes_dict[s_index] = ("cone", p_index)

            elif s.type == "torus":
                extplanes, inplanes = surf_planes
                p_ext = []
                p_in = []
                for p in extplanes:
                    planes[next_index] = p
                    p_ext.append(next_index)
                    next_index += 1
                for p in inplanes:
                    planes[next_index] = p
                    p_in.append(next_index)
                    next_index += 1
                surf_planes_dict[s_index] = ("torus", p_ext, p_in)
            else:
                p_index = []
                for p in surf_planes:
                    planes[next_index] = p
                    p_index.append(next_index)
                    next_index += 1
                surf_planes_dict[s_index] = p_index
    return planes, surf_planes_dict


def convert_to_planes(s, pos):
    if s.type == "cylinder":
        return cylinder_to_planes(s, pos)
    elif s.type == "cone":
        return cone_to_planes(s, pos)
    elif s.type == "sphere":
        return sphere_to_planes(s, pos)
    elif s.type == "torus":
        return torus_to_planes(s, pos)
    elif s.type == "paraboloid":
        return parabola_to_planes(s, pos)
    else:
        print(f"{s.type} not implemented for boundbox")
        return []


def get_orto_axis(axis):
    x = FreeCAD.Vector(1, 0, 0)
    z = FreeCAD.Vector(0, 0, 1)
    vx = axis.cross(x)
    vz = axis.cross(z)
    if vx.Length < vz.Length:
        v = vz
    else:
        v = vx
    v.normalize()
    w = v.cross(axis)
    w.normalize()

    return v, w


def cylinder_to_planes(cyl, pos):
    center, axis, radius = cyl.params
    if pos is None:
        radius = radius * 0.8535533906
    elif pos:
        radius = radius * 0.70710678

    x, y = get_orto_axis(axis)
    r1 = center + x * radius
    r2 = center - x * radius
    r3 = center + y * radius
    r4 = center - y * radius

    p1 = Part.Plane(r1, -x)
    p2 = Part.Plane(r2, x)
    p3 = Part.Plane(r3, -y)
    p4 = Part.Plane(r4, y)
    return (p1, p2, p3, p4)


def cone_to_planes(cone, pos):
    apex, axis, t, dbl = cone.params
    if pos is None:
        t = t * 0.8535533906
    elif pos:
        t = t * 0.70710678
    sa = math.atan(t)
    nface = 4
    x, y = get_orto_axis(axis)
    cs = math.cos(sa)
    ss = math.sin(sa)
    dphi = twoPi / nface
    phi = 0
    cplanes = []
    for i in range(nface):
        rho = x * math.cos(phi) + y * math.sin(phi)
        ni = -axis * ss + rho * cs
        pi = Part.Plane(apex, -ni)
        cplanes.append(pi)
        phi += dphi

    pa = Part.Plane(apex, axis)
    cplanes.append(pa)

    return cplanes


def sphere_to_planes(sphere, pos):
    center, radius = sphere.params
    if pos is None:
        radius = radius * 0.8535533906
    elif pos:
        radius = radius * 0.70710678

    x = FreeCAD.Vector(1, 0, 0)
    y = FreeCAD.Vector(0, 1, 0)
    z = FreeCAD.Vector(0, 0, 1)

    r1 = center + x * radius
    r2 = center - x * radius
    r3 = center + y * radius
    r4 = center - y * radius
    r5 = center + z * radius
    r6 = center - z * radius

    p1 = Part.Plane(r1, -x)
    p2 = Part.Plane(r2, x)
    p3 = Part.Plane(r3, -y)
    p4 = Part.Plane(r4, y)
    p5 = Part.Plane(r5, -z)
    p6 = Part.Plane(r6, z)
    return (p1, p2, p3, p4, p5, p6)


def torus_to_planes(torus, pos):
    center, axis, majorRadius, minorR, minorA = torus.params

    x = FreeCAD.Vector(1, 0, 0)
    y = FreeCAD.Vector(0, 1, 0)
    z = FreeCAD.Vector(0, 0, 1)

    if pos is None:
        dist = (majorRadius + minorR) * 0.8535533906
        difR = (majorRadius - minorR) * 0.8535533906
    elif pos:
        dist = (majorRadius + minorR) * 0.70710678
        difR = majorRadius - minorR
    else:
        dist = majorRadius + minorR
        difR = (majorRadius - minorR) * 0.70710678

    if abs(abs(axis.dot(x)) - 1) < 1e-5:
        r1 = center + x * minorA
        r2 = center - x * minorA
        r3 = center + y * dist
        r4 = center - y * dist
        r5 = center + z * dist
        r6 = center - z * dist
        if difR > 0:
            r7 = center + y * difR
            r8 = center - y * difR
            r9 = center + z * difR
            r10 = center - z * difR
            p7 = Part.Plane(r7, y)
            p8 = Part.Plane(r8, -y)
            p9 = Part.Plane(r9, z)
            p10 = Part.Plane(r10, -z)
    elif abs(abs(axis.dot(y)) - 1) < 1e-5:
        r1 = center + x * dist
        r2 = center - x * dist
        r3 = center + y * minorA
        r4 = center - y * minorA
        r5 = center + z * dist
        r6 = center - z * dist
        if difR > 0:
            r7 = center + x * difR
            r8 = center - x * difR
            r9 = center + z * difR
            r10 = center - z * difR
            p7 = Part.Plane(r7, x)
            p8 = Part.Plane(r8, -x)
            p9 = Part.Plane(r9, z)
            p10 = Part.Plane(r10, -z)
    elif abs(abs(axis.dot(z)) - 1) < 1e-5:
        r1 = center + x * dist
        r2 = center - x * dist
        r3 = center + y * dist
        r4 = center - y * dist
        r5 = center + z * minorA
        r6 = center - z * minorA
        if difR > 0:
            r7 = center + x * difR
            r8 = center - x * difR
            r9 = center + y * difR
            r10 = center - y * difR
            p7 = Part.Plane(r7, x)
            p8 = Part.Plane(r8, -x)
            p9 = Part.Plane(r9, y)
            p10 = Part.Plane(r10, -y)

    p1 = Part.Plane(r1, -x)
    p2 = Part.Plane(r2, x)
    p3 = Part.Plane(r3, -y)
    p4 = Part.Plane(r4, y)
    p5 = Part.Plane(r5, -z)
    p6 = Part.Plane(r6, z)
    external_planes = (p1, p2, p3, p4, p5, p6)
    if difR > 0:
        central_planes = (p7, p8, p9, p10)
    else:
        central_planes = tuple()
    return (external_planes, central_planes)


def parabola_to_planes(parabola, pos):
    # parabola approximated by plane tanget to the curve
    # plane separation such that distance from plane to curve < a*x0 (a parameter < 1, x0 absica of tangent point )
    # the sequence of tangent points is xn+1 = xn * (1 - sqrt(2*a))
    # initial point x0 is calculated suche that after n iteration the last y = xn^2 / (4*focal) is < rmax, where rmax is the maximin universe distance
    # x0 = b^n * sqrt(4*focal*rmax)  whith b =  (1 - sqrt(2*a); n number of tangent planes to consider
    #

    center, axis, focal = parabola.params
    nt = 7
    a = 0.22
    rmax = 7.5e5
    b = 1 - math.sqrt(2 * a)
    x0 = b**nt * math.sqrt(4 * focal * rmax)

    axis.normalize()
    x, y = get_orto_axis(axis)
    dphi = twoPi / 4
    phi = 0

    p0 = Part.Plane(center, axis)
    cplanes = [p0]
    focal = float(focal)
    xp = x0
    for n in range(nt + 1):
        zi = 0.25 * xp * xp / focal
        for i in range(4):
            rho = x * math.cos(phi) + y * math.sin(phi)
            vec = y * math.cos(phi) - x * math.sin(phi)
            slope = 2 * focal * rho + xp * axis  # slope xp/(2*focal)
            xe = center + xp * rho + zi * axis
            normal = vec.cross(slope)
            normal.normalize()
            pi = Part.Plane(xe, normal)
            cplanes.append(pi)
            phi += dphi
        xp = xp / b

    return cplanes


def plane_definition(seq, surf_index, orientation):
    for s, planes in surf_index.items():
        if len(planes) == 0:
            continue
        addpos = False
        if type(planes[0]) is str:
            if planes[0] == "torus":
                extplanes, inplanes = planes[1:3]
                extm = BoolSequence(" ".join((str(p) for p in extplanes)))
                if len(inplanes) > 0:
                    inm = BoolSequence(":".join((str(p) for p in inplanes)))
                    pm = BoolSequence(operator="AND")
                    pm.append(extm, inm)
                else:
                    pm = extm
            elif planes[0] == "dblcone":
                addpos = True
                cplanes = planes[1]
                cone1 = BoolSequence(" ".join((str(p) for p in cplanes)))
                cone2 = BoolSequence(" ".join((str(-p) for p in cplanes)))
                pm = BoolSequence(operator="OR")
                pm.append(cone1, cone2)
            else:
                cplanes = planes[1]
                addpos = True
                pm = BoolSequence(" ".join((str(p) for p in cplanes)))
        else:
            pm = BoolSequence(" ".join((str(p) for p in planes)))
        if type(seq.elements) is bool:
            return seq
        pp = pm.get_complementary()
        if orientation == "Forward":
            change_surf(seq, -s, pm)
            if addpos:
                change_surf(seq, s, pp)
            else:
                change_surf(seq, s, True)
        elif orientation == "Reversed":
            change_surf(seq, s, pp)
            if addpos:
                change_surf(seq, -s, pm)
            else:
                change_surf(seq, -s, False)
        else:
            change_surf(seq, s, pp)
            change_surf(seq, -s, pm)

    seq.join_operators()
    return seq


def change_surf(seq, old, new):
    clean = False
    for i, e in enumerate(seq.elements):
        if type(e) is BoolSequence:
            if abs(old) in e.get_surfaces_numbers():
                change_surf(e, old, new)
                if type(e.elements) is bool:
                    if e.elements == (seq.operator == "OR"):
                        seq.elements = e.elements
                        return
                    else:
                        clean = True
        else:
            if e == old:
                if type(new) is bool:
                    if new == (seq.operator == "OR"):
                        seq.elements = new
                        return
                    else:
                        seq.elements[i] = new
                        clean = True
                else:
                    seq.elements[i] = new
    if clean:
        seq.clean()


def plane_intersect(plane_list, externalBox, cutBoundary):
    point_list = []
    if not cutBoundary:
        for i, p1 in enumerate(plane_list[0:-2]):
            j = i + 1
            for p2 in plane_list[i + 1 : -1]:
                line = p1.intersect(p2)
                if len(line) == 0:
                    continue
                line = line[0]
                for p3 in plane_list[j + 1 :]:
                    inter = line.intersect(p3)
                    if len(inter[0]) == 0:
                        continue
                    p = inter[0][0]
                    p = FreeCAD.Vector(p.X, p.Y, p.Z)
                    if externalBox.isInside(p):
                        point_list.append(p)
                j += 1
    else:
        XYZ = (
            FreeCAD.Vector(1, 0, 0),
            FreeCAD.Vector(0, 1, 0),
            FreeCAD.Vector(0, 0, 1),
        )
        pxm = Part.Plane(XYZ[0], FreeCAD.Vector(externalBox.XMin, 0, 0))
        pxp = Part.Plane(XYZ[0], FreeCAD.Vector(externalBox.XMax, 0, 0))
        pym = Part.Plane(XYZ[1], FreeCAD.Vector(0, externalBox.YMin, 0))
        pyp = Part.Plane(XYZ[1], FreeCAD.Vector(0, externalBox.YMax, 0))
        pzm = Part.Plane(XYZ[2], FreeCAD.Vector(0, 0, externalBox.ZMin))
        pzp = Part.Plane(XYZ[2], FreeCAD.Vector(0, 0, externalBox.ZMax))
        PXYZ = (pxm, pxp, pym, pyp, pzm, pzp)

        for i, p1 in enumerate(plane_list[0:]):
            j = i + 1
            point_list.extend(plane_boundary(p1, externalBox))
            for p2 in plane_list[i + 1 :]:
                line = p1.intersect(p2)
                if len(line) == 0:
                    continue
                line = line[0]
                point_list.extend(line_boundary(line, externalBox, PXYZ))
                for p3 in plane_list[j + 1 :]:
                    inter = line.intersect(p3)
                    if len(inter[0]) == 0:
                        continue
                    p = inter[0][0]
                    p = FreeCAD.Vector(p.X, p.Y, p.Z)
                    if externalBox.isInside(p):
                        point_list.append(p)

        for i in range(8):
            p = externalBox.getPoint(i)
            point_list.append(p)

    return point_list


def plane_boundary(plane, externalBox):

    point_list = []
    for i in range(12):
        segment = Part.LineSegment(*externalBox.getEdge(i))
        inter = segment.intersect(plane)
        if len(inter[0]) == 1:
            p = inter[0][0]
            p = FreeCAD.Vector(p.X, p.Y, p.Z)
            point_list.append(p)
    return point_list


def line_boundary(line, externalBox, PXYZ):
    points = []
    for i, plane in enumerate(PXYZ):
        inter = line.intersect(plane)
        if len(inter[0]) == 0:
            continue
        p = inter[0][0]
        p = FreeCAD.Vector(p.X, p.Y, p.Z)
        if i < 2:
            if (externalBox.YMin <= p.y <= externalBox.YMax) and (externalBox.ZMin <= p.z <= externalBox.ZMax):
                points.append(p)
        elif i < 4:
            if (externalBox.XMin <= p.x <= externalBox.XMax) and (externalBox.ZMin <= p.z <= externalBox.ZMax):
                points.append(p)
        else:
            if (externalBox.XMin <= p.x <= externalBox.XMax) and (externalBox.YMin <= p.y <= externalBox.YMax):
                points.append(p)
        if len(points) == 2:
            return points
    return points


def sort_point(point_list, axis):
    axis_points = []
    if axis == "x":
        for i, point in enumerate(point_list):
            axis_points.append((point.x, i))
    elif axis == "y":
        for i, point in enumerate(point_list):
            axis_points.append((point.y, i))
    elif axis == "z":
        for i, point in enumerate(point_list):
            axis_points.append((point.z, i))
    else:
        print("bad axis name")

    axis_points.sort()
    sorted_points = list(point_list[x[1]] for x in axis_points)
    removed = remove_close_points(list(sorted_points))
    return removed


def remove_close_points(sorted_list):
    if len(sorted_list) < 2:
        return sorted_list
    new_points = []
    p = sorted_list.pop()
    new_points.append(p)
    while len(sorted_list) > 0:
        nextp = sorted_list.pop()
        dp = p - nextp
        while dp.Length < 0.1:
            if len(sorted_list) > 0:
                nextp = sorted_list.pop()
                dp = p - nextp
            else:
                break
        else:
            p = nextp
            new_points.append(p)
    new_points.reverse()
    return new_points


def remove_points(point_list, value, axis, lower, Lmax=None):
    kept_points = []
    if axis == "x":
        if lower:
            for p in point_list[::-1]:
                if p.x < value:
                    break
                kept_points.append(p)
        else:
            for p in point_list[::-1]:
                if p.x > value:
                    break
                kept_points.append(p)
    elif axis == "y":
        if lower:
            for p in point_list[::-1]:
                if p.y < value:
                    break
                kept_points.append(p)
        else:
            for p in point_list[::-1]:
                if p.y > value:
                    break
                kept_points.append(p)
    elif axis == "z":
        if lower:
            for p in point_list[::-1]:
                if p.z < value:
                    break
                kept_points.append(p)
        else:
            for p in point_list[::-1]:
                if p.z > value:
                    break
                kept_points.append(p)
    else:
        print("bad axis name")
    return kept_points


def pointaxis(p, axis):
    if axis == "x":
        return p.x
    elif axis == "y":
        return p.y
    elif axis == "z":
        return p.z


def makePlane(normal, position, Box):

    p0 = normal.dot(position)

    pointEdge = []
    for i in range(12):
        edge = Box.getEdge(i)
        p1 = normal.dot(edge[0])
        p2 = normal.dot(edge[1])
        d0 = p0 - p1
        d1 = p2 - p1
        if d1 != 0:
            a = d0 / d1
            if a >= 0 and a <= 1:
                pointEdge.append(edge[0] + a * (edge[1] - edge[0]))

    if len(pointEdge) == 0:
        return None  # Plane does not cross box

    s = FreeCAD.Vector((0, 0, 0))
    for v in pointEdge:
        s = s + v
    s = s / len(pointEdge)

    vtxvec = []
    for v in pointEdge:
        vtxvec.append(v - s)

    X0 = vtxvec[0]
    Y0 = normal.cross(X0)

    orden = []
    for i, v in enumerate(vtxvec):
        phi = numpy.arctan2(v.dot(Y0), v.dot(X0))
        orden.append((phi, i))
    orden.sort()

    return Part.Face(Part.makePolygon([pointEdge[p[1]] for p in orden], True))


def inertia_matrix(points):
    npoints = len(points)
    numpy_points = numpy.ndarray((npoints, 3))
    for i, p in enumerate(points):
        numpy_points[i] = numpy.array((p.x, p.y, p.z))

    x0 = numpy.sum(numpy_points[:, 0]) / npoints
    y0 = numpy.sum(numpy_points[:, 1]) / npoints
    z0 = numpy.sum(numpy_points[:, 2]) / npoints
    Sxx = numpy.sum(numpy_points[:, 0] * numpy_points[:, 0])
    Syy = numpy.sum(numpy_points[:, 1] * numpy_points[:, 1])
    Szz = numpy.sum(numpy_points[:, 2] * numpy_points[:, 2])
    Sxy = numpy.sum(numpy_points[:, 0] * numpy_points[:, 1])
    Sxz = numpy.sum(numpy_points[:, 0] * numpy_points[:, 2])
    Syz = numpy.sum(numpy_points[:, 1] * numpy_points[:, 2])

    Ixx = Sxx / npoints - x0 * x0
    Iyy = Syy / npoints - y0 * y0
    Izz = Szz / npoints - z0 * z0
    Ixy = Sxy / npoints - x0 * y0
    Ixz = Sxz / npoints - x0 * z0
    Iyz = Syz / npoints - y0 * z0

    inertia_matrix = numpy.array(((Ixx, Ixy, Ixz), (Ixy, Iyy, Iyz), (Ixz, Iyz, Izz)))
    eigvalue, vectors = numpy.linalg.eig(inertia_matrix)
    return


def box_intersect(Fbox, Rbox):
    PX1 = (Fbox.Box.XMin, Fbox.Box.XMax)
    PX2 = (Rbox.Box.XMin, Rbox.Box.XMax)
    PY1 = (Fbox.Box.YMin, Fbox.Box.YMax)
    PY2 = (Rbox.Box.YMin, Rbox.Box.YMax)
    PZ1 = (Fbox.Box.ZMin, Fbox.Box.ZMax)
    PZ2 = (Rbox.Box.ZMin, Rbox.Box.ZMax)

    orientation = Fbox.Orientation
    bXmin, bXmax = Fbox.Box.XMin, Fbox.Box.XMax
    bYmin, bYmax = Fbox.Box.YMin, Fbox.Box.YMax
    bZmin, bZmax = Fbox.Box.ZMin, Fbox.Box.ZMax

    xmin, xmax = plane_region(PX1, PX2, orientation)
    boxes = []
    if xmin is not None:
        box = FreeCAD.BoundBox(xmin, bYmin, bZmin, xmax, bYmax, bZmax)
        boxes.append(box)

    ymin, ymax = plane_region(PY1, PY2, orientation)
    if ymin is not None:
        box = FreeCAD.BoundBox(bXmin, ymin, bZmin, bXmax, ymax, bZmax)
        boxes.append(box)

    zmin, zmax = plane_region(PZ1, PZ2, orientation)
    if zmin is not None:
        box = FreeCAD.BoundBox(bXmin, bYmin, zmin, bXmax, bYmax, zmax)
        boxes.append(box)

    if len(boxes) > 0:
        box = boxes[0]
        for b in boxes[1:]:
            box.add(b)
        return box
    else:
        return None


def plane_region(P1, P2, orient1):
    p11, p12 = P1
    p21, p22 = P2

    if p11 >= p22:
        return (p11, p12) if orient1 == "Forward" else (p21, p22)
    elif p12 <= p21:
        return (p11, p12) if orient1 == "Forward" else (p21, p22)
    else:
        if p11 < p21:
            if p12 < p22:
                return (p11, p21) if orient1 == "Forward" else (p12, p22)
            else:
                return (p11, p12) if orient1 == "Forward" else (None, None)
        elif p11 > p21:
            if p12 <= p22:
                return (None, None) if orient1 == "Forward" else (p21, p22)  # OK
            else:
                return (p22, p12) if orient1 == "Forward" else (p21, p11)
        else:
            if p12 < p22:
                return (None, None) if orient1 == "Forward" else (p12, p22)
            elif p12 > p22:
                return (p22, p12) if orient1 == "Forward" else (None, None)
            else:
                return (None, None)


def operate_box(definition, boxes):
    fullbox = myBox()
    definition.level_update()
    for e in definition.elements:
        if type(e) is int:
            box = boxes[abs(e)]
            if definition.operator == "AND":
                fullbox.mult(box)
            else:
                fullbox.add(box)
        else:
            box = operate_box(e, boxes)
            if definition.operator == "AND":
                fullbox.mult(box)
            else:
                fullbox.add(box)
    return fullbox
