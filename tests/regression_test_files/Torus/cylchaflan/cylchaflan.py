# openMC geometry script generated by GEOUNED
import openmc

###############################################################################
# Define problem geometry
###############################################################################


# Surface setup
S1 = openmc.YPlane(y0=0.0)
S2 = openmc.ZPlane(z0=0.0)
S3 = openmc.ZPlane(z0=10.0)
S4 = openmc.ZPlane(z0=10.200000000000001)
S5 = openmc.ZCylinder(x0=0.0,y0=0.0,r=2.0)
S6 = openmc.ZTorus(x0=0.0,y0=0.0,z0=10.0,r=1.8,r1=0.2,r2=0.2)
S7 = openmc.ZCone(x0=0.0,y0=0.0,z0=8.4,r2=0.9999999999999998)
S8 = openmc.XPlane(x0=-3.164784400584788)
S9 = openmc.XPlane(x0=3.164784400584788)
S10 = openmc.YPlane(y0=-3.164784400584788)
S11 = openmc.YPlane(y0=3.164784400584788)
S12 = openmc.ZPlane(z0=-1.0)
S13 = openmc.ZPlane(z0=11.200000000000001)
S14 = openmc.Sphere(x0=0.0,y0=0.0,z0=5.1000000000000005,r=7.717142354316536, boundary_type="vacuum")

# Cell definition 
C1 = openmc.Cell(name="", region=((+S2 & -S5 & -S3) | (-S6 & +S3) | (-S7 & +S3 & -S4 & +S6)))
C2 = openmc.Cell(name="Automatic Generated Void Cell. Enclosure(-3.165, 3.165, -3.165, 3.165, -1.000, 11.200). Enclosed cells : (1)", region=(+S8 & -S9 & +S10 & -S11 & +S12 & -S13 & ((+S6 | -S3) & (+S5 | -S2 | +S3) & (+S7 | +S4 | -S3 | -S6))))
C3 = openmc.Cell(name="", region=(-S14 & (-S8 | +S9 | -S10 | +S11 | -S12 | +S13)))

geometry = openmc.Geometry([C1, C2, C3])
geometry.export_to_xml()
