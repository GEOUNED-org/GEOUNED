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
        verbose (bool, optional): Print output warning during geouned run.
            Defaults to False.
        enlargeBox (float, optional): Enlarge box boundary when evaluating
            the constraint table during the simplification of the void cell
            definition. (unit is millimeter). Defaults to 2.
        nPlaneReverse (int, optional): Threshold value to determine whether
            cut with parallel planes should be carried out first. Defaults
            to 0.
        splitTolerance (float, optional): Fuzzy tolerance value used in the
            FreeCAD function “BOPTools.SplitAPI.slice”. This function is
            used during the solid decomposition process. Defaults to 0.
        scale_up (bool, optional): Scale up Fuzzy tolerance once get below
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
        verbose: bool = False,
        enlargeBox: float = 2.0,
        nPlaneReverse: int = 0,
        splitTolerance: float = 0.0,
        scale_up: bool = True,
        quadricPY: bool = False,
        Facets: bool = False,
        prnt3PPlane: bool = False,
        forceNoOverlap: bool = False,
    ):

        self.forceCylinder = forceCylinder
        self.newSplitPlane = newSplitPlane
        self.delLastNumber = delLastNumber
        self.verbose = verbose
        self.enlargeBox = enlargeBox
        self.nPlaneReverse = nPlaneReverse
        self.splitTolerance = splitTolerance
        self.scale_up = scale_up
        self.quadricPY = quadricPY
        self.Facets = Facets
        self.prnt3PPlane = prnt3PPlane
        self.forceNoOverlap = forceNoOverlap


class Tolerances:
    """_summary_

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
        min_area (float, optional): minimun face area to consider in cell definition. Defaults to 1.0e-2.
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


class NumericFormat:
    """Numerical format options for each of the surface types.

    Args:
        P_abc (str, optional): Plane general a,b,c params. Defaults to "14.7e".
        P_d (str, optional): Plane general d params. Defaults to "14.7e".
        P_xyz (str, optional): PX/PY/PZ params. Defaults to "14.7e".
        S_r (str, optional): SO/SX/SY/SZ/S radius. Defaults to "14.7e".
        S_xyz (str, optional): SO/SX/SY/SZ/S center. Defaults to "14.7e".
        C_r (str, optional): Cylinder radius. Defaults to "12f".
        C_xyz (str, optional): Cylinder center. Defaults to "12f".
        K_xyz (str, optional): Cone apex. Defaults to "13.6e".
        K_tan2 (str, optional): Cone tan^2 value. Defaults to "12f".
        T_r (str, optional): Torus radii. Defaults to "14.7e".
        T_xyz (str, optional): Torus center. Defaults to "14.7e".
        GQ_1to6 (str, optional): GQ 1 to 6 coefficients (order 2 x2,y2,z2,xy,...). Defaults to "18.15f".
        GQ_7to9 (str, optional): GQ 7 to 9 coefficients (order 1 x,y,z). Defaults to "18.15f".
        GQ_10 (str, optional): GQ 10 coefficient. Defaults to "18.15f".
    """

    def __init__(
        self,
        P_abc: str = "14.7e",
        P_d: str = "14.7e",
        P_xyz: str = "14.7e",
        S_r: str = "14.7e",
        S_xyz: str = "14.7e",
        C_r: str = "12f",
        C_xyz: str = "12f",
        K_xyz: str = "13.6e",
        K_tan2: str = "12f",
        T_r: str = "14.7e",
        T_xyz: str = "14.7e",
        GQ_1to6: str = "18.15f",
        GQ_7to9: str = "18.15f",
        GQ_10: str = "18.15f",
    ):

        self.P_abc = P_abc
        self.P_d = P_d
        self.P_xyz = P_xyz
        self.S_r = S_r
        self.S_xyz = S_xyz
        self.C_r = C_r
        self.C_xyz = C_xyz
        self.K_xyz = K_xyz
        self.K_tan2 = K_tan2
        self.T_r = T_r
        self.T_xyz = T_xyz
        self.GQ_1to6 = GQ_1to6
        self.GQ_7to9 = GQ_7to9
        self.GQ_10 = GQ_10


class Settings:
    """Settings for changing the way the CAD to CSG conversion is done

    Args:
        stepFile (str, optional): Name of the CAD file (in STEP format) to
            be converted. Defaults to "".
        matFile (str, optional): _description_. Defaults to "".
        voidGen (bool, optional): Generate voids of the geometry. Defaults
            to True.
        debug (bool, optional): Write step files of original and decomposed
            solids, for each solid in the STEP file. Defaults to False.
        compSolids (bool, optional): Join subsolids of STEP file as a single
            compound solid. Step files generated with SpaceClaim have not
            exactly the same level of solids as FreeCAD. It may a happened
            that solids defined has separated solids are read by FreeCAD
            as a single compound solid (and will produce only one MCNP
            cell). In this case compSolids should be set to False. Defaults
            to True.
        simplify (str, optional): Simplify the cell definition considering
            relative surfaces position and using Boolean logics. Available
            options are: "no" no optimization, "void" only void cells are
            simplified. Algorithm is faster but the simplification is not
            optimal. "voidfull" : only void cells are simplified with the
            most optimal algorithm. The time of the conversion can be
            multiplied by 5 or more. "full" : all the cells (solids and
            voids) are simplified. Defaults to "No".
        cellRange (list, optional): Range of cell to be converted (only one
            range is allowed, e.g [100,220]). Default all solids are
            converted. Defaults to [].
        exportSolids (str, optional): Export CAD solid after reading.
            The execution is stopped after export, the translation is not
            carried out. Defaults to "".
        minVoidSize (float, optional): Minimum size of the edges of the
            void cell. Units are in mm. Defaults to 200.0.
        max_surf (int, optional): #TODO
        max_bracket (int, optional): Maximum number of brackets (solid
            complementary) allowed in void cell definition. Defaults to 30.
        voidMat (list, optional): Assign a material defined by the user
            instead of void for cells without material definition and the
            cells generated in the automatic void generation. The format
            is a 3 valued tuple (mat_label, mat_density, mat_description).
            Example (100,1e-3,'Air assigned to Void'). Defaults to [].
        voidExclude (list, optional): #TODO see issue 87. Defaults to [].
        startCell (int, optional): Starting cell numbering label. Defaults to 1.
        startSurf (int, optional): Starting surface numbering label. Defaults to 1.
        sort_enclosure (bool, optional): If enclosures are defined in the
            CAD models, the voids cells of the enclosure will be located in
            the output file in the same location where the enclosure solid
            is located in the CAD solid tree.. Defaults to False.
    """

    def __init__(
        self,
        matFile: str = "",
        voidGen: bool = True,
        debug: bool = False,
        compSolids: bool = True,
        simplify: str = "No",
        cellRange=[],
        exportSolids: str = "",
        minVoidSize: float = 200.0,  # units mm
        max_surf: int = 50,
        max_bracket: int = 30,
        voidMat=[],
        voidExclude=[],
        startCell: int = 1,
        startSurf: int = 1,
        sort_enclosure: bool = False,
    ):

        self.matFile = matFile
        self.voidGen = voidGen
        self.debug = debug
        self.compSolids = compSolids
        self.simplify = simplify
        self.cellRange = cellRange
        self.exportSolids = exportSolids
        self.minVoidSize = minVoidSize
        self.max_surf = max_surf
        self.max_bracket = max_bracket
        self.voidMat = voidMat
        self.voidExclude = voidExclude
        self.startCell = startCell
        self.startSurf = startSurf
        self.sort_enclosure = sort_enclosure