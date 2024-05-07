from pathlib import Path


import geouned

path_to_cad = Path("testing/inputSTEP")
step_files = list(path_to_cad.rglob("*.stp")) + list(path_to_cad.rglob("*.step"))

input_step_file = step_files[0]

# sets up an output folder for the result

my_options = geouned.Options(
    forceCylinder=True,
    splitTolerance=0,
    newSplitPlane=True,
    nPlaneReverse=0,
)

my_settings = geouned.Settings(
    compSolids=False,
    voidGen=True,
    minVoidSize=100,
    debug=False,
    simplify="no",
)

geo = geouned.CadToCsg(
    step_file=f"{input_step_file.resolve()}", settings=my_settings, options=my_options
)

geo.start()

geo.export_csg(
    title="Input Test",
    geometry_name="test",
    out_formats=("mcnp", "openmc_xml", "openmc_py", "serpent", "phits"),
    volCARD=False,
    volSDEF=True,
    cellSummaryFile=False,
    cellCommentFile=False,
    dummyMat=True,
)
