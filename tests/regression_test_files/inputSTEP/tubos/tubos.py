# openMC geometry script generated by GEOUNED
import openmc

###############################################################################
# Define problem geometry
###############################################################################


# Surface setup
S1 = openmc.Plane(a=0.4689449996294938,b=-0.4282426204033891,c=0.7724628440206901,d=-1289.638833198076)
S2 = openmc.Plane(a=-0.5294241686781376,b=0.8423604593834029,c=-0.100692135178699,d=1852.4828400982199)
S3 = openmc.Cylinder(x0=-961.7978028047636,y0=1491.541250906148,z0=-326.16443172000004,r=3.49965,dx=0.7244847831937575,dy=-0.6616034336585792,dz=0.19339776496090208)
S4 = openmc.Plane(a=0.2923717047223038,b=-0.956304755963168,c=8.159775107656208e-13,d=-1795.0805769740155)
S5 = openmc.Cylinder(x0=-1009.2941440671559,y0=1550.8027954178106,z0=-337.044472998,r=3.49965,dx=0.2923717047223038,dy=-0.956304755963168,dz=8.163417574676206e-13)
S6 = openmc.ZPlane(z0=-277.9998680728499)
S7 = openmc.ZCylinder(x0=-912.9498129378808,y0=1446.9330217072306,r=3.49965)
S8 = openmc.XPlane(x0=-1018.5969970333416)
S9 = openmc.XPlane(x0=-907.4530344047115)
S10 = openmc.YPlane(y0=1441.4083363198881)
S11 = openmc.YPlane(y0=1569.0367354915766)
S12 = openmc.ZPlane(z0=-341.6693716211191)
S13 = openmc.ZPlane(z0=-276.9999979988159)
S14 = openmc.Sphere(x0=-963.0250157190266,y0=1505.2225359057325,z0=-309.3346848099675,r=92.39887837360851, boundary_type="vacuum")

# Cell definition 
C1 = openmc.Cell(name="", region=((-S3 & -S1 & -S2)))
C2 = openmc.Cell(name="", region=((+S2 & -S5 & +S4)))
C3 = openmc.Cell(name="", region=((-S7 & -S6 & +S1)))
C4 = openmc.Cell(name="Automatic Generated Void Cell. Enclosure(-1018.597, -907.453, 1441.408, 1569.037, -341.669, -277.000). Enclosed cells : (1, 2, 3)", region=(+S8 & -S9 & +S10 & -S11 & +S12 & -S13 & ((+S3 | +S1 | +S2) & (+S5 | -S4 | -S2) & (+S7 | +S6 | -S1))))
C5 = openmc.Cell(name="", region=(-S14 & (-S8 | +S9 | -S10 | +S11 | -S12 | +S13)))

geometry = openmc.Geometry([C1, C2, C3, C4, C5])
geometry.export_to_xml()