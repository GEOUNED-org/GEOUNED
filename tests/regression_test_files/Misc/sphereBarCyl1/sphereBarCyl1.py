# openMC geometry script generated by GEOUNED
import openmc

###############################################################################
# Define problem geometry
###############################################################################


# Surface setup
S1 = openmc.YPlane(y0=0.30000000000000004)
S2 = openmc.Plane(a=0.7071067811865476,b=-0.0,c=0.7071067811865476,d=1.2121320343558972)
S3 = openmc.Cylinder(x0=0.1,y0=0.2,z0=0.2,r=0.2,dx=0.7071067811865476,dy=0.0,dz=0.7071067811865476)
S4 = openmc.Sphere(x0=0.0,y0=0.30000000000000004,z0=0.1,r=0.5)
S5 = openmc.Plane(a=0.6899402905873934,b=-0.2190086598684918,c=0.689940288957523,d=0.0032914309352047645)
S6 = openmc.XPlane(x0=-0.8806078918683564)
S7 = openmc.XPlane(x0=1.9485281374238597)
S8 = openmc.YPlane(y0=-1.0696367414507255)
S9 = openmc.YPlane(y0=1.4324128622454548)
S10 = openmc.ZPlane(z0=-0.7586291902410247)
S11 = openmc.ZPlane(z0=2.0485281374238595)
S12 = openmc.Sphere(x0=0.5339601227777516,y0=0.18138806039736466,z0=0.6449494735914175,r=2.3999494118665328, boundary_type="vacuum")

# Cell definition 
C1 = openmc.Cell(name="", region=((+S5 & +S4 & -S3 & -S2)))
C2 = openmc.Cell(name="Automatic Generated Void Cell. Enclosure(-0.881, 1.949, -1.070, 1.432, -0.759, 2.049). Enclosed cells : (1)", region=(+S6 & -S7 & +S8 & -S9 & +S10 & -S11 & (-S5 | -S4 | +S3 | +S2)))
C3 = openmc.Cell(name="", region=(-S12 & (-S6 | +S7 | -S8 | +S9 | -S10 | +S11)))

geometry = openmc.Geometry([C1, C2, C3])
geometry.export_to_xml()