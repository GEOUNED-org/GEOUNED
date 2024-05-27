# openMC geometry script generated by GEOUNED
import openmc

###############################################################################
# Define problem geometry
###############################################################################


# Surface setup
S1 = openmc.Plane(a=-0.6030848728438236,b=0.10634013473996888,c=0.7905570263367687,d=-701.9643762043311)
S2 = openmc.Cylinder(x0=720.8296916985,y0=-102.22377115100001,z0=-324.8476130149,r=12.5,dx=0.9010177524985242,dy=-0.1588737398629161,dz=-0.40364110849178686)
S3 = openmc.Cylinder(x0=765.282773405934,y0=-110.06204883170652,z0=-344.76186123744355,r=9.999999999999602,dx=-0.9010177524985242,dy=0.1588737398629161,dz=0.40364110849178686)
S4 = openmc.Cylinder(x0=718.8048357018001,y0=-101.8667344068,z0=-323.94051088870003,r=7.499999999999701,dx=-0.9010177524985242,dy=0.1588737398629161,dz=0.40364110849178686)
S5 = openmc.ZCone(x0=0.0,y0=0.0,z0=-1971.2852599682953,r2=0.2093159503611188)
S6 = openmc.ZCone(x0=0.0,y0=0.0,z0=-2004.9361661889359,r2=0.2093159503611188)
S7 = openmc.XPlane(x0=707.8605085936322)
S8 = openmc.XPlane(x0=760.6074892340785)
S9 = openmc.YPlane(y0=-120.05661731457135)
S10 = openmc.YPlane(y0=-87.72813872822752)
S11 = openmc.ZPlane(z0=-350.02554602804213)
S12 = openmc.ZPlane(z0=-309.47823970257195)
S13 = openmc.Sphere(x0=734.2339989138553,y0=-103.89237802139945,z0=-329.75189286530707,r=37.7243471442727, boundary_type="vacuum")

# Cell definition 
C1 = openmc.Cell(name="", region=((+S3 & -S5 & -S2 & -S1) | (-S6 & +S4 & -S3 & -S1)))
C2 = openmc.Cell(name="Automatic Generated Void Cell. Enclosure(707.861, 760.607, -120.057, -87.728, -350.026, -309.478). Enclosed cells : (1)", region=(+S7 & -S8 & +S9 & -S10 & +S11 & -S12 & ((+S1) | ((+S5 | +S2 | -S3) & (+S6 | -S4 | +S3)))))
C3 = openmc.Cell(name="", region=(-S13 & (-S7 | +S8 | -S9 | +S10 | -S11 | +S12)))

geometry = openmc.Geometry([C1, C2, C3])
geometry.export_to_xml()