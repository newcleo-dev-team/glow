"""
Testing the construction of a hexagonal cell. The 'MATERIAL' property is
assigned to each region of the cell technological geometry.
"""
import os
import sys

from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.geometries import Circle
from glow.geometry_layouts.layouts import Region
from glow.support.types import LayoutGeometryType, PropertyType
from glow.main import TdtSetup, export_layout_to_tdt


edge_length = 0.7852193995
# Build a hexagonal cell filled with 'LEAD' and with four circular regions
hex_cell = HexCell(
    side=edge_length,
    name="Hexagonal cell",
    base_props={PropertyType.MATERIAL: "LEAD"}
)
for r, mat in zip(
        [0.525, 0.465, 0.45, 0.1],
        ["CLADDING", "HOLLOW", "FUEL", "HOLLOW"]
    ):
    hex_cell.add(
        Region(Circle(radius=r), properties={PropertyType.MATERIAL: mat})
    )

# Display the cell's regions with the 'MATERIAL' property type colour map
hex_cell.show(PropertyType.MATERIAL)

# Generate the output TDT file from the cell
export_layout_to_tdt(
    hex_cell,
    os.path.join(os.path.dirname(sys.argv[0]), 'test_hex_cell'),
    TdtSetup(type_geo=LayoutGeometryType.HEXAGON_TRAN)
)
