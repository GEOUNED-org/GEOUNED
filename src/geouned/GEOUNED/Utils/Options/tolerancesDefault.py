defaultValues = {
    "relativeTol": False,
    "relativePrecision": 1.0e-6,  # relative precision
    "value": 1.0e-6,  # Tolerance in single value comparison
    "distance": 1.0e-4,  # General Distance Tolerance
    "angle": 1.0e-4,  # General Angle Tolerance
    "pln_distance": 1.0e-4,  # distance between planes equal planes if distance between parallel planes < 1e-4 cm
    "pln_angle": 1.0e-4,  # angle between axis. 1e-4 : planes separate each other 0.1mm each 1m
    "cyl_distance": 1.0e-4,  # distance between radius/center
    "cyl_angle": 1.0e-4,  # angle between axis
    "sph_distance": 1.0e-4,  # distance between radius/center
    "kne_distance": 1.0e-4,  # distance between apex
    "kne_angle": 1.0e-4,  # angle between semiangles/axis
    "tor_distance": 1.0e-4,  # distance between Major/Minor radii/center
    "tor_angle": 1.0e-4,  # angle between axis
    "min_area": 1.0e-2,  # minimun face area to consider in cell definition
}

typeDict = {
    "relativeTol": bool,
    "relativePrecision": float,  # relative precision
    "value": float,  # Tolerance in single value comparison
    "distance": float,  # General Distance Tolerance
    "angle": float,  # General Angle Tolerance
    "pln_distance": float,  # distance between planes equal planes if distance between parallel planes < 1e-4 cm
    "pln_angle": float,  # angle between axis. 1e-4 : planes separate each other 0.1mm each 1m
    "cyl_distance": float,  # distance between radius/center
    "cyl_angle": float,  # angle between axis
    "sph_distance": float,  # distance between radius/center
    "kne_distance": float,  # distance between apex
    "kne_angle": float,  # angle between semiangles/axis
    "tor_distance": float,  # distance between Major/Minor radii/center
    "tor_angle": float,  # angle between axis
    "min_area": float,  # minimun face area to consider in cell definition
}

KwrdEquiv = {
    "relativeTolerance": "relativeTol",
    "relativePrecision": "relativePrecision",
    "singleValue": "value",
    "generalDistance": "distance",
    "generalAngle": "angle",
    "planeDistance": "pln_distance",
    "planeAngle": "pln_angle",
    "cylinderDistance": "cyl_distance",
    "cylinderAngle": "cyl_angle",
    "sphereDistance": "sph_distance",
    "coneDistance": "kne_distance",
    "coneAngle": "kne_angle",
    "torusDistance": "tor_distance",
    "torusAngle": "tor_angle",
    "minArea": "min_area",
}
