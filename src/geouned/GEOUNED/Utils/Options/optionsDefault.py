# default options values
defaultValues = {
    "forceCylinder": False,  # Force using cylinders instead cones for auxillary surfaces of torus surface definition
    "newSplitPlane": True,  # Use new module for planes splitting in decomposition module
    "delLastNumber": False,  # Deleting the last word in the comment if it is a number
    "verbose": False,  # Display information during conversion process
    "enlargeBox": 2,  # Enlarge box Boundary when evaluating Ctable (dimension in mm)
    "nPlaneReverse": 0,  # numbers of plane thresold whether cut of parallel planes are made fisrt with lowest or highest number
    "splitTolerance": 0,  # First try for fuzzy tolerance used in BOPTOOL Split function. If BOPTSplit crash try lower value for tolerance
    "scaleUp": True,  # Scale up Fuzzy tolerance once get below 1e-12
    "quadricPY": False,  # use quadric form of cones and cylinder not aligned with X,Y,Z axis when write openMC script file
    "Facets": False,  # use alternative conversion module when geometry is defined by cells compound by only triangular plane faces
    "prnt3PPlane": False,  # print 3 point plane definition in MCNP output as 3 points coordinates
    "forceNoOverlap": False,  # force no overlaping cell definition. Adjacent cell definition are rested from current cell definition
}

typeDict = {
    "forceCylinder": bool,  # Force using cylinders instead cones for auxillary surfaces of torus surface definition
    "newSplitPlane": bool,  # Use new module for planes splitting in decomposition module
    "delLastNumber": bool,  # Deleting the last word in the comment if it is a number
    "verbose": bool,  # Display information during conversion process
    "enlargeBox": float,  # Enlarge box Boundary when evaluating Ctable (dimension in mm)
    "nPlaneReverse": int,  # numbers of plane thresold whether cut of parallel planes are made fisrt with lowest or highest number
    "splitTolerance": int,  # First try for fuzzy tolerance used in BOPTOOL Split function. If BOPTSplit crash try lower value for tolerance
    "scaleUp": bool,  # Scale up Fuzzy tolerance once get below 1e-12
    "quadricPY": bool,  # use quadric form of cones and cylinder not aligned with X,Y,Z axis when write openMC script file
    "Facets": bool,  # use alternative conversion module when geometry is defined by cells compound by only triangular plane faces
    "prnt3PPlane": bool,  # print 3 point plane definition in MCNP output as 3 points coordinates
    "forceNoOverlap": bool,  # force no overlaping cell definition. Adjacent cell definition are rested from current cell definition
}
