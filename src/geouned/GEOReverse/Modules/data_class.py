class Options:
    splitTolerance = 1.0e-2


class BoxSettings:
    """Parameters used in the solids boundbox generation. Optimized dimensions can reduce
    the translation time.

    Args:
        universe_radius (float, optional): Maximum radius of the CAD universe.
            Solids with coordinates x^2+y^2+z*2 > universe_radius^2 will be cut or not represented.
            Units mm. Defaults to 1.0e6.
        insolid_tolerance (float, optional): Maximum distance from the nearest
            surface of the solid, for which a point outside the solid is assumed
            inside the solid. Used only for boundbox generation. Units mm.
            Defaults to 1.
    """

    def __init__(
        self,
        universe_radius: float = 1.0e6,  # units mm
        insolid_tolerance: float = 1,  # units mm
    ):

        self.universe_radius = universe_radius
        self.insolid_tolerance = insolid_tolerance

    @property
    def universe_radius(self):
        return self._universe_radius

    @universe_radius.setter
    def universe_radius(self, universe_radius: float):
        if not isinstance(universe_radius, (float, int)):
            raise TypeError(f"geoReverse.Settings.universe_radius should be a float, not a {type(universe_radius)}")
        self._universe_radius = universe_radius

    @property
    def insolid_tolerance(self):
        return self._insolid_tolerance

    @insolid_tolerance.setter
    def insolid_tolerance(self, insolid_tolerance: float):
        if not isinstance(insolid_tolerance, (float, int)):
            raise TypeError(f"geoReverse.Settings.insolid_tolerance should be a float, not a {type(insolid_tolerance)}")
        self._insolid_tolerance = insolid_tolerance
