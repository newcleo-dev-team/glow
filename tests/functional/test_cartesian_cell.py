"""
Testing the construction of a cartesian cell. The 'MATERIAL' property is
assigned to each region of the cell technological geometry.
"""
import os
import sys

from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.geometries import Circle
from glow.geometry_layouts.layouts import Region
from glow.support.types import LayoutGeometryType, PropertyType
from glow.main import TdtSetup, export_layout_to_tdt


# Build a Cartesian cell filled with 'COOLANT' and with three circular regions
rect_cell = CartesianCell(
    name="Cartesian cell",
    base_props={PropertyType.MATERIAL: "COOLANT"}
)
# Add three inner circles to the cell
radii = [0.65/2, 0.3, 0.75/2]
for r, mat in zip([0.375, 0.325, 0.3], ["CLADDING", "HOLLOW", "FUEL"]):
    rect_cell.add(
        Region(Circle(radius=r), properties={PropertyType.MATERIAL: mat})
    )

# Display the cell's regions with the 'MATERIAL' property type colour map
rect_cell.show(PropertyType.MATERIAL)

# Generate the output TDT file from the cell
export_layout_to_tdt(
    rect_cell,
    os.path.join(os.path.dirname(sys.argv[0]), 'test_cartesian_cell'),
    TdtSetup(type_geo=LayoutGeometryType.RECTANGLE_TRAN)
)