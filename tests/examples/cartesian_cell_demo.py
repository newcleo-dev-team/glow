"""
Testing the construction of a cartesian cell. The 'MATERIAL' property is
assigned to each region of the cell technological geometry. A sectorization
is applied to the cell and the result graphically shown in the SALOME 3D
viewer.
"""
from glow.geometry_layouts.cells import RectCell
from glow.generator.support import GeometryType, PropertyType
from glow.geometry_layouts.geometries import Circle
from glow.interface.geom_interface import ShapeType, get_point_coordinates, \
    make_partition, make_vertex

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
c1 = Circle(center=get_point_coordinates(v1), radius=0.075)
c2 = Circle(center=get_point_coordinates(v2), radius=0.075)
c3 = Circle(center=get_point_coordinates(v3), radius=0.075)
c4 = Circle(center=get_point_coordinates(v4), radius=0.075)
c1.build_face()
c2.build_face()
c3.build_face()
c4.build_face()
# Partition the original cell face with the 4 circles
updated_face = make_partition(
    [rect_cell.face, c1.face, c2.face, c3.face, c4.face], [], ShapeType.FACE)
# Update the cell geometric layout with the just built shape
rect_cell.update_geometry_from_face(GeometryType.TECHNOLOGICAL, updated_face)
# Show the result in the 3D viewer
rect_cell.show(PropertyType.MATERIAL)
