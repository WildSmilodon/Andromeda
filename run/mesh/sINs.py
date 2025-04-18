import gmsh

gmsh.initialize()
gmsh.model.add("sphere_in_sphere")

# Parameters
sizeSmallSphere = 0.1
radiusSmallSphere = 0.5
sizeLageSphere = 150.0
radiusLargeSphere = 512.0

# Small sphere
points_small = [
    gmsh.model.geo.addPoint(0, 0, 0, sizeSmallSphere),
    gmsh.model.geo.addPoint(-radiusSmallSphere, 0, 0, sizeSmallSphere),
    gmsh.model.geo.addPoint(radiusSmallSphere, 0, 0, sizeSmallSphere),
    gmsh.model.geo.addPoint(0, -radiusSmallSphere, 0, sizeSmallSphere),
    gmsh.model.geo.addPoint(0, radiusSmallSphere, 0, sizeSmallSphere),
    gmsh.model.geo.addPoint(0, 0, radiusSmallSphere, sizeSmallSphere),
    gmsh.model.geo.addPoint(0, 0, -radiusSmallSphere, sizeSmallSphere)
]

circles_small = [
    gmsh.model.geo.addCircleArc(points_small[2], points_small[0], points_small[4]),
    gmsh.model.geo.addCircleArc(points_small[1], points_small[0], points_small[4]),
    gmsh.model.geo.addCircleArc(points_small[1], points_small[0], points_small[3]),
    gmsh.model.geo.addCircleArc(points_small[3], points_small[0], points_small[2]),
    gmsh.model.geo.addCircleArc(points_small[5], points_small[0], points_small[2]),
    gmsh.model.geo.addCircleArc(points_small[5], points_small[0], points_small[1]),
    gmsh.model.geo.addCircleArc(points_small[1], points_small[0], points_small[6]),
    gmsh.model.geo.addCircleArc(points_small[6], points_small[0], points_small[2]),
    gmsh.model.geo.addCircleArc(points_small[4], points_small[0], points_small[5]),
    gmsh.model.geo.addCircleArc(points_small[5], points_small[0], points_small[3]),
    gmsh.model.geo.addCircleArc(points_small[3], points_small[0], points_small[6]),
    gmsh.model.geo.addCircleArc(points_small[6], points_small[0], points_small[4])
]

line_loops_small = [
    gmsh.model.geo.addCurveLoop([circles_small[8], circles_small[4], circles_small[0]]),
    gmsh.model.geo.addCurveLoop([circles_small[1], circles_small[8], circles_small[5]]),
    gmsh.model.geo.addCurveLoop([circles_small[6], circles_small[11], -circles_small[1]]),
    gmsh.model.geo.addCurveLoop([circles_small[11], -circles_small[0], -circles_small[7]]),
    gmsh.model.geo.addCurveLoop([circles_small[5], circles_small[2], -circles_small[9]]),
    gmsh.model.geo.addCurveLoop([circles_small[9], circles_small[3], -circles_small[4]]),
    gmsh.model.geo.addCurveLoop([circles_small[10], -circles_small[6], circles_small[2]]),
    gmsh.model.geo.addCurveLoop([circles_small[7], -circles_small[3], circles_small[10]])
]

surfaces_small = [gmsh.model.geo.addSurfaceFilling([loop]) for loop in line_loops_small]


gmsh.model.geo.synchronize()  # Synchronize before adding physical groups
gmsh.model.addPhysicalGroup(2, surfaces_small, 1)
gmsh.model.setPhysicalName(2, 1, "particle")

## Large sphere
points_large = [
    gmsh.model.geo.addPoint(0, 0, 0, sizeLageSphere),
    gmsh.model.geo.addPoint(-radiusLargeSphere, 0, 0, sizeLageSphere),
    gmsh.model.geo.addPoint(radiusLargeSphere, 0, 0, sizeLageSphere),
    gmsh.model.geo.addPoint(0, -radiusLargeSphere, 0, sizeLageSphere),
    gmsh.model.geo.addPoint(0, radiusLargeSphere, 0, sizeLageSphere),
    gmsh.model.geo.addPoint(0, 0, radiusLargeSphere, sizeLageSphere),
    gmsh.model.geo.addPoint(0, 0, -radiusLargeSphere, sizeLageSphere)
]

circles_large = [
    gmsh.model.geo.addCircleArc(points_large[2], points_large[0], points_large[4]),
    gmsh.model.geo.addCircleArc(points_large[1], points_large[0], points_large[4]),
    gmsh.model.geo.addCircleArc(points_large[1], points_large[0], points_large[3]),
    gmsh.model.geo.addCircleArc(points_large[3], points_large[0], points_large[2]),
    gmsh.model.geo.addCircleArc(points_large[5], points_large[0], points_large[2]),
    gmsh.model.geo.addCircleArc(points_large[5], points_large[0], points_large[1]),
    gmsh.model.geo.addCircleArc(points_large[1], points_large[0], points_large[6]),
    gmsh.model.geo.addCircleArc(points_large[6], points_large[0], points_large[2]),
    gmsh.model.geo.addCircleArc(points_large[4], points_large[0], points_large[5]),
    gmsh.model.geo.addCircleArc(points_large[5], points_large[0], points_large[3]),
    gmsh.model.geo.addCircleArc(points_large[3], points_large[0], points_large[6]),
    gmsh.model.geo.addCircleArc(points_large[6], points_large[0], points_large[4])
]

line_loops_large = [
    gmsh.model.geo.addCurveLoop([circles_large[8], circles_large[4], circles_large[0]]),
    gmsh.model.geo.addCurveLoop([circles_large[1], circles_large[8], circles_large[5]]),
    gmsh.model.geo.addCurveLoop([circles_large[6], circles_large[11], -circles_large[1]]),
    gmsh.model.geo.addCurveLoop([circles_large[11], -circles_large[0], -circles_large[7]]),
    gmsh.model.geo.addCurveLoop([circles_large[5], circles_large[2], -circles_large[9]]),
    gmsh.model.geo.addCurveLoop([circles_large[9], circles_large[3], -circles_large[4]]),
    gmsh.model.geo.addCurveLoop([circles_large[10], -circles_large[6], circles_large[2]]),
    gmsh.model.geo.addCurveLoop([circles_large[7], -circles_large[3], circles_large[10]])
]

surfaces_large = [gmsh.model.geo.addSurfaceFilling([loop]) for loop in line_loops_large]


gmsh.model.geo.synchronize()  # Synchronize before adding physical groups
gmsh.model.addPhysicalGroup(2, surfaces_large, 2)
gmsh.model.setPhysicalName(2, 2, "outside")


gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)  # Generate mesh first

# Reverse normals for specific surfaces
for surface in [surfaces_large[1], surfaces_large[2], surfaces_large[6]]:
    gmsh.model.mesh.reverse([(2, surface)])  # Reverse normals after meshing

# Reverse normals for specific surfaces
for surface in [surfaces_small[4], surfaces_small[3], surfaces_small[0], surfaces_small[5], surfaces_small[7]]:
    gmsh.model.mesh.reverse([(2, surface)])  # Reverse normals after meshing


gmsh.write("sINs-2D.msh")  # Gmsh format
gmsh.write("sINs-2D.vtk")  # VTK format


## Define mesh size fields
field_small = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(field_small, "NodesList", points_small)
#
field_large = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(field_large, "NodesList", points_large)
#
threshold_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(threshold_field, "InField", field_small)
gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", sizeSmallSphere)
gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", sizeLageSphere)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", radiusSmallSphere)
gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", radiusLargeSphere)
gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)


## Create surface loops and volume
surface_loop_large = gmsh.model.geo.addSurfaceLoop(surfaces_large)
surface_loop_small = gmsh.model.geo.addSurfaceLoop(surfaces_small)
volume = gmsh.model.geo.addVolume([surface_loop_large, surface_loop_small])



gmsh.model.geo.synchronize()  # Synchronize before adding physical groups
gmsh.model.addPhysicalGroup(3, [volume], 3)
gmsh.model.setPhysicalName(3, 3, "fluid")

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("sINs-3D.msh")  # Gmsh format
gmsh.write("sINs-3D.vtk")  # VTK format

# Write mesh to file
# Save mesh in different formats
#gmsh.write("sINS.msh")  # Gmsh format
#gmsh.write("sINS.vtk")  # VTK format
#gmsh.write("sINS.stl")  # STL format
#gmsh.write("sINS.unv")  # UNV format
#gmsh.write("sINS.med")  # MED format
#gmsh.write("sINS.cgns")  # CGNS format
#gmsh.write("sINS.ply")  # PLY format
#gmsh.write("sINS.obj")  # OBJ format

gmsh.finalize()