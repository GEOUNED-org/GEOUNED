class Options:
    """A class for containing conversion options

    Args:
        forceCylinder (bool, optional): Use cylinder (instead of cones) as
            ancillary surface where unclosed torus surfaces are involved in
            the solid definition. Defaults to False.
        newSplitPlane (bool, optional): New method to consider plane as
            cutting surface during the decomposition process. Former method
            split first planes perpendicular to X,Y,Z axis and then the
            other planes involved in the solid definition. New method group
            all parallel planes independently whether their normal are
            along X,Y,Z axes, and start the decomposition process cutting
            first with the group having the highest number of parallel
            planes. Defaults to True.
        delLastNumber (bool, optional): Deleting the last word in the
            comment if it is a number. Defaults to False.
        enlargeBox (float, optional): Enlarge box boundary when evaluating
            the constraint table during the simplification of the void cell
            definition. (unit is millimeter). Defaults to 2.
        nPlaneReverse (int, optional): Threshold value to determine whether
            cut with parallel planes should be carried out first. Defaults
            to 0.
        splitTolerance (float, optional): Fuzzy tolerance value used in the
            FreeCAD function “BOPTools.SplitAPI.slice”. This function is
            used during the solid decomposition process. Defaults to 0.
        scaleUp (bool, optional): Scale up Fuzzy tolerance once get below
            1e-12. Defaults to True.
        quadricPY (bool, optional): In openMC python script format, the
            cones or cylinders no aligned with the X,Y, or Z axis can be
            defined using the openmc.Cone or open.Cylinder methods but can
            also be defined with their quadric parameter. If “quadricPY” is
            11 True then all cones and cylinders will be defined in the
            openMC python script format under their quadric form. Defaults
            to False.
        Facets (bool, optional): use alternative conversion module when
            geometry is defined by cells compound by only triangular plane
            faces. Defaults to False.
        prnt3PPlane (bool, optional): print 3 point plane definition in
            output as 3 points coordinates. Defaults to False.
        forceNoOverlap (bool, optional): force no overlaping cell
            definition. Adjacent cell definition are rested from current
            cell definition. Defaults to False.
    """

    def __init__(
        self,
        forceCylinder: bool = False,
        newSplitPlane: bool = True,
        delLastNumber: bool = False,
        enlargeBox: float = 2.0,
        nPlaneReverse: int = 0,
        splitTolerance: float = 0.0,
        scaleUp: bool = True,
        quadricPY: bool = False,
        Facets: bool = False,
        prnt3PPlane: bool = False,
        forceNoOverlap: bool = False,
    ):

        self.forceCylinder = forceCylinder
        self.newSplitPlane = newSplitPlane
        self.delLastNumber = delLastNumber
        self.enlargeBox = enlargeBox
        self.nPlaneReverse = nPlaneReverse
        self.splitTolerance = splitTolerance
        self.scaleUp = scaleUp
        self.quadricPY = quadricPY
        self.Facets = Facets
        self.prnt3PPlane = prnt3PPlane
        self.forceNoOverlap = forceNoOverlap


class Tolerances:
    """A class for containing tolerances values

    Args:
        relativeTol (bool, optional): _description_. Defaults to False.
        relativePrecision (float, optional): relative precision. Defaults to 1.0e-6.
        value (float, optional): Tolerance in single value comparison. Defaults to 1.0e-6.
        distance (float, optional): General Distance Tolerance. Defaults to 1.0e-4.
        angle (float, optional): General Angle Tolerance. Defaults to 1.0e-4.
        pln_distance (float, optional): distance between planes equal planes if distance between parallel planes < 1e-4 cm. Defaults to 1.0e-4.
        pln_angle (float, optional): angle between axis. 1e-4 : planes separate each other 0.1mm each 1m. Defaults to 1.0e-4.
        cyl_distance (float, optional): distance between radius/center. Defaults to 1.0e-4.
        cyl_angle (float, optional): angle between axis. Defaults to 1.0e-4.
        sph_distance (float, optional): distance between radius/center. Defaults to 1.0e-4.
        kne_distance (float, optional): distance between apex. Defaults to 1.0e-4.
        kne_angle (float, optional): angle between semiangles/axis. Defaults to 1.0e-4.
        tor_distance (float, optional): distance between Major/Minor radii/center. Defaults to 1.0e-4.
        tor_angle (float, optional): angle between axis. Defaults to 1.0e-4.
        min_area (float, optional): minimum face area to consider in cell definition. Defaults to 1.0e-2.
    """

    def __init__(
        self,
        relativeTol: bool = False,
        relativePrecision: float = 1.0e-6,
        value: float = 1.0e-6,
        distance: float = 1.0e-4,
        angle: float = 1.0e-4,
        pln_distance: float = 1.0e-4,
        pln_angle: float = 1.0e-4,
        cyl_distance: float = 1.0e-4,
        cyl_angle: float = 1.0e-4,
        sph_distance: float = 1.0e-4,
        kne_distance: float = 1.0e-4,
        kne_angle: float = 1.0e-4,
        tor_distance: float = 1.0e-4,
        tor_angle: float = 1.0e-4,
        min_area: float = 1.0e-2,
    ):

        self.relativeTol = relativeTol
        self.relativePrecision = relativePrecision
        self.value = value
        self.distance = distance
        self.angle = angle
        self.pln_distance = pln_distance
        self.pln_angle = pln_angle
        self.cyl_distance = cyl_distance
        self.cyl_angle = cyl_angle
        self.sph_distance = sph_distance
        self.kne_distance = kne_distance
        self.kne_angle = kne_angle
        self.tor_distance = tor_distance
        self.tor_angle = tor_angle
        self.min_area = min_area
