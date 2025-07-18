"""
Testing the construction of a cartesian cell. The 'MATERIAL' property is
assigned to each region of the cell technological geometry. A sectorization
is applied to the cell and the result graphically shown in the SALOME 3D
viewer.
"""
from glow.geometry_layouts.cells import RectCell
from glow.support.types import GeometryType, PropertyType
from glow.geometry_layouts.geometries import Circle
from glow.interface.geom_interface import *

# -----------------------
# TEST THE CARTESIAN CELL
# -----------------------
# Build a cartesian cell
rect_cell = RectCell(name="Cartesian cell")
# Add three inner circles to the cell
radii = [0.65/2, 0.3, 0.75/2]
for r in radii:
    rect_cell.add_circle(r)

# Assign the materials to each zone in the cell
rect_cell.set_properties(
    {PropertyType.MATERIAL: ["TSTR.'TMIL_MOC'.'COMB0201'.1",
                             "TSTR.'TMIL_MOC'.'COMB0201'.2",
                             "TSTR.'TMIL_MOC'.'COMB0201'.3",
                             "TSTR.'TMIL_MOC'.'MODEXTCRAYON'"]}
)

# Build the cell sectorization
rect_cell.sectorize([1, 1, 4, 8], [0, 0, 0, 22.5], windmill=True)
# Show the sectorized cell with regions colored according to the 'MATERIAL'
# property
rect_cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)
# Print the properties assigned to each cell sector
for sec in rect_cell.regions:
    print(f"Sector region {sec.name}, {sec.properties}")

# -----------------------------------------------------------------------
# Test the possibility to update the cell face from a given face built in
# SALOME.
# -----------------------------------------------------------------------
# Build vertex objects representing the centers of the circles to add
v1 = make_vertex((0, 0.15, 0))
v2 = make_vertex((0, -0.15, 0))
v3 = make_vertex((0.15, 0, 0))
v4 = make_vertex((-0.15, 0, 0))
# Build the corresponding 'Circle' objects
circles = [Circle(center=get_point_coordinates(v),
                  radius=0.075) for v in [v1, v2, v3, v4]]
# Build circles positioned in the cell center
center_circles = [Circle(radius=r) for r in [0.1, 0.35, 0.45]]
# Update the list of 'Circle' objects
circles += center_circles

# Partition the original cell face with all the circles
updated_face = make_partition(
    [rect_cell.face], [c.face for c in circles], ShapeType.FACE)
# Update the cell geometric layout with the just built shape
rect_cell.update_geometry_from_face(GeometryType.TECHNOLOGICAL, updated_face)
# Show the result in the 3D viewer
rect_cell.show(PropertyType.MATERIAL)
