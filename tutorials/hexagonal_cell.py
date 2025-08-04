"""
Use case showing the construction of a hexagonal cell whose technological
geometry is made by four regions.
The 'MATERIAL' property is assigned to each region of the cell's technological
geometry. A sectorization is applied to the cell and the result graphically
shown in the SALOME 3D viewer.
"""
from glow.geometry_layouts.cells import HexCell
from glow.support.types import GeometryType, PropertyType


# Build the cell's geometry layout by adding three circular regions
cell = HexCell(name="Hexagonal Cell")
radii = [0.25, 0.4, 0.6]
for radius in radii:
    cell.add_circle(radius)
# Show the cell's technological geometry in the SALOME 3D viewer
cell.show()

# Assign the materials to each zone in the cell
cell.set_properties(
      {PropertyType.MATERIAL: ["MAT_1", "MAT_2", "MAT_3", "MAT_4"]}
)
# Show the regions by applying a colorset
cell.show(PropertyType.MATERIAL)

# Build the cell's sectorized geometry
cell.sectorize([1, 1, 6, 6], [0]*4)
# Show the sectorized cell with regions colored according to the 'MATERIAL'
# property
cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)
