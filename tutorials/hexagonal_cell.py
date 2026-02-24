"""
Use case showing the construction of a hexagonal cell whose technological
geometry is made by four regions.
The 'MATERIAL' property is assigned to each region of the cell's technological
geometry. A sectorization is applied to the cell and the result graphically
shown in the SALOME 3D viewer.
"""
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.geometries import Circle
from glow.geometry_layouts.layouts import Region
from glow.support.types import GeometryType, PropertyType


# Intialise three lists, one storing the circular regions radii, the other the
# names of the materials, sorted from the inner to the outer region
radii = [0.25, 0.4, 0.6]
materials = ["MAT_1", "MAT_2", "MAT_3"]
# Build the cell's geometry layout by adding three circular regions from the
# outer to the inner
cell = HexCell(
    name="Cartesian cell", base_props={PropertyType.MATERIAL: "MAT_4"}
)
for radius, mat in zip(radii[::-1], materials[::-1]):
    cell.add(
        Region(Circle(radius=radius), properties={PropertyType.MATERIAL: mat})
    )

# Show the regions according to the 'MATERIAL' colour map
cell.show(PropertyType.MATERIAL)

# Build the cell's sectorized geometry
cell.sectorize([1, 1, 6, 6], [0]*4)
# Show the cell's sectorized layout according to the 'MATERIAL' colour map
cell.show(PropertyType.MATERIAL, GeometryType.SECTORIZED)
