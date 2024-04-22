from .Classes import Options, Tolerances, MCNP_numeric_format

#default options values
forceCylinder  = True                           # Force using cylinders instead cones for auxillary surfaces of torus surface definition
newSplitPlane  = True                           # Use new module for planes splitting in decomposition module
delLastNumber  = False                          # Deleting the last word in the comment if it is a number
verbose        = False                          # Display information during conversion process
enlargeBox     = 2                              # Enlarge box Boundary when evaluating Ctable (dimension in mm)
nPlaneReverse  = 0                              # numbers of plane thresold whether cut of parallel planes are made fisrt with lowest or highest number
splitTolerance = 0                              # First try for fuzzy tolerance used in BOPTOOL Split function. If BOPTSplit crash try lower value for tolerance
scaleUp        = True                           # Scale up Fuzzy tolerance once get below 1e-12
quadricPY      = False                          # use quadric form of cones and cylinder not aligned with X,Y,Z axis when write openMC script file
Facets         = False                          # use alternative conversion module when geometry is defined by cells compound by only triangular plane faces
prnt3PPlane    = False                          # print 3 point plane definition in MCNP output as 3 points coordinates
forceNoOverlap = False                          # force no overlaping cell definition. Adjacent cell definition are rested from current cell definition

tolValueDict = { 'relativeTol'        : False ,      
              'relativePrecision'  : 1.e-6 ,   # relative precision
              'value'              : 1.e-6 ,   # Tolerance in single value comparison
              'distance'           : 1.e-4 ,   # General Distance Tolerance
              'angle'              : 1.e-4 ,   # General Angle Tolerance
              'pln_distance'       : 1.e-4 ,   # distance between planes equal planes if distance between parallel planes < 1e-4 cm
              'pln_angle'          : 1.e-4 ,   # angle between axis. 1e-4 : planes separate each other 0.1mm each 1m 
              'cyl_distance'       : 1.e-4 ,   # distance between radius/center
              'cyl_angle'          : 1.e-4 ,   # angle between axis
              'sph_distance'       : 1.e-4 ,   # distance between radius/center
              'kne_distance'       : 1.e-4 ,   # distance between apex
              'kne_angle'          : 1.e-4 ,   # angle between semiangles/axis
              'tor_distance'       : 1.e-4 ,   # distance between Major/Minor radii/center
              'tor_angle'          : 1.e-4 ,   # angle between axis
              'min_area'           : 1.e-2 }   # minimun face area to consider in cell definition


numValueDict = {  'P_abc'   : '14.7e'  ,     # Plane general a,b,c params
                  'P_d'     : '14.7e'  ,     # Plane general d     params
                  'P_xyz'   : '14.7e'  ,     # PX/PY/PZ            params
                  'S_r'     : '14.7e'  ,     # SO/SX/SY/SZ/S       radius
                  'S_xyz'   : '14.7e'  ,     # SO/SX/SY/SZ/S       center 
                  'C_r'     : '12f'    ,     # Cylinder            radius
                  'C_xyz'   : '12f'    ,     # Cylinder            center 
                  'K_xyz'   : '13.6e'  ,     # Cone                apex
                  'K_tan2'  : '12f'    ,     # Cone                tan^2 value
                  'T_r'     : '14.7e'  ,     # Torus               radii
                  'T_xyz'   : '14.7e'  ,     # Torus               center 
                  'GQ_1to6' : '18.15f' ,     # GQ 1 to 6 coefficients (order 2 x2,y2,z2,xy,...)
                  'GQ_7to9' : '18.15f' ,     # GQ 7 to 9 coefficients (order 1 x,y,z)
                  'GQ_10'   : '18.15f' }     # GQ 10 coefficient      (order 0)



# Set default  attributes to Options class
setattr(Options,'forceCylinder' ,forceCylinder)  
setattr(Options,'newSplitPlane' ,newSplitPlane)  
setattr(Options,'enlargeBox'    ,enlargeBox)  
setattr(Options,'delLastNumber' ,delLastNumber)  
setattr(Options,'verbose'       ,verbose)  
setattr(Options,'nPlaneReverse' ,nPlaneReverse)  
setattr(Options,'splitTolerance',splitTolerance)  
setattr(Options,'scaleUp'       ,scaleUp)  
setattr(Options,'quadricPY'     ,quadricPY)  
setattr(Options,'Facets'        ,Facets)  
setattr(Options,'prnt3PPlane'   ,prnt3PPlane)  
setattr(Options,'forceNoOverlap',forceNoOverlap)  

# Set default  attributes to Tolerances class
for key in tolValueDict.keys():
    setattr(Tolerances,key,tolValueDict[key])


# Set default  attributes to MCNP number format class
for key in numValueDict.keys():
    setattr(MCNP_numeric_format,key,numValueDict[key])
