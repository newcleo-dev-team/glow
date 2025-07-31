"""
Testing the construction of a hexagonal cell. The 'MATERIAL' property is
assigned to each region of the cell technological geometry. A sectorization
is applied to the cell and the result graphically shown in the SALOME 3D
viewer.
"""
# -----------------------
# TEST THE HEXAGONAL CELL
# -----------------------
from glow.geometry_layouts.cells import HexCell
from glow.support.types import GeometryType, PropertyType


hex_cell = HexCell(name="Hexagonal cell")
# Rotate the cell
radii = [0.65/2, 0.3, 0.75/2]
# Add three inner circles to the cell
for r in radii:
    hex_cell.add_circle(r)

# Add both cell face and its subfaces from the technological geometry
# to the SALOME study
hex_cell.show()

# Assign the materials to each zone in the cell
hex_cell.set_properties(
    {PropertyType.MATERIAL: ["TSTR.'TMIL_MOC'.'COMB0201'.1",
                             "TSTR.'TMIL_MOC'.'COMB0201'.2",
                             "TSTR.'TMIL_MOC'.'COMB0201'.3",
                             "TSTR.'TMIL_MOC'.'MODEXTCRAYON'"]}
)

# Build the cell sectorization
hex_cell.sectorize([1, 1, 6, 6], [0]*4)

# Show the sectorized cell with regions colored according to the 'MATERIAL'
# property
hex_cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)
# Print the properties assigned to each cell sector
for sec in hex_cell.regions:
    print(f"Sector region {sec.name}, {sec.properties}")
