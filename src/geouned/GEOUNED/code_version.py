""" Modification from 0.9.5 
  - volInput card name substituted by volSDEF card name
  - new card volCARD (values : True/False) used to insert the CAD cell volume into the cell definition
    using the MCNP card VOL. Default value True  
  - Fixed bug. Code crashing when automatic void calculation was not enabled.
  - Fixed bug. Comment line might not be in the correct order following de CAD tree labels

Modification in 0.9.6 Release Date 24/11/2022
  - Small change in the reading of material file.
  - In "CellDefinition.py" file. "is_opposite" function is call with rhe correct tolerance value instead of default value of the is_opposite function.  

Modification in 0.9.7 Release Date 06/12/2022
  - Compatibility with FreeCAD version higher than 0.18.
  - fixed bug with sphere Objects

Modification in 0.9.7 Release Date 06/12/2022
  - Removed external packages Boolean. Subtituted by own make package
  - Full cell optimization using CTable approach and own factorization function
  - Cylinders instead cones can be used  as auxillary surface for Torus surface definition 
  - fix bugs several additional planes for the same reversed surface
  - new module for plane splitting in decomposition modules. All planes are considerer equivalent (PX,PY,PZ,P), parallel planes are grouped together. 
    Decomposition is performed by splitting first with group of parallel planes with lowest number of elements. 
    Old version module is still available by switching the option keyword "newSplitPlane" to False. 
  - New option added when enclosure are used. The option "sort_enclosure" True (Default False) will produce output file where solids and void of an enclosure
    are grouped together (first solid cell of the enclosure, then voids)
  - Options has been added to enlarge Box dimension when Constraint table is evaluated
 
 Modification in 0.9.8 
  - fix bug with two cylinders/cones intersection
  - Introduce splitting planes to cut rounded edges(edge made with cylinder)
  - fixed bug with torus (I hope so!)
  - new card "dummyMat" to write dummy material card (associated to material labels present in cell definition) in input file.

 new in version 0.9.8.2
  - Javi modification for memory bug for geomtry with large plane numbers
  - small changes in decomposition process
  - new verbose option to remove some output on screen

 new in version 1.0
  - allow to convert geometry to openMC XML format and script input file
  - change the entry "mcnpfile" by "geometryName" in FILES section. Now geometryName is the generic name of output geometry file.
    mncp input file will have ".mcnp" extension and openMC XML the ".xml" extension.
  - new entry "outFormat" is added in FILES section. "outFormat" is used to select the output geometry format (mcnp and/or openMC_xml). 
    Default value is only mcnp
  - new entry splitTolerance. Change the Tolerance value of the FreeCAD function split. Default Value  0.
  - new entry Verbose. Print on screen warning. Default False.

 new in version 1.0.1
  - implement new conversion module when all cells are made of triangular plane faces
  - bug fixed void complementary of level1 'AND' sequence
  - changes in simplify function of BoolSequence class
  - implementation of the geometry conversion to Serpent code input file (by A. Valentine UKAEA)
  - fixed bug with torodial surfaces.
  - implementation of new, more robust, cell definition, where adjacent solid cell definition is rested to current solid cell definition (forcing no overlaping cell).
  - fix bug additional plane with sphere
  - improvement of cell simplification method
  - fix bug. bad parameter normalization multiple output formats
  - Change tolerance values if BOPTsplit function crash
  - fix some bugs with Torus (not all)
  - implementation of the geometry conversion to PHITS code input file (by T. Okano Nissin High Voltage Corporation (NHVC))
"""
