"""
Testing the construction of a cartesian cell. The 'MATERIAL' property is
assigned to each region of the cell technological geometry. A sectorization
is applied to the cell and the result graphically shown in the SALOME 3D
viewer.
"""
from glow.geometry_layouts.cells import RectCell
from glow.generator.support import GeometryType, PropertyType

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
