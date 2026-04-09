"""
Testing the construction of a lattice made by cartesian cells. The 'MATERIAL'
property is assigned to each region of the cell technological geometry.
A lattice is built by including 36 of the built cells. Symmetry is exploited,
one-eighth of the complete cartesian lattice is considered.
"""
import os
import sys

from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.geometries import Circle
from glow.geometry_layouts.layouts import Region
from glow.support.types import *
from glow.geometry_layouts.lattices import CartesianLattice
from glow.main import TdtSetup, export_layout_to_tdt


# Build a Cartesian cell filled with 'COOLANT' and with three circular regions
rect_cell = CartesianCell(
    name="Cartesian cell",
    base_props={PropertyType.MATERIAL: "COOLANT"}
)
# Add three inner circles to the cell
for r, mat in zip([0.375, 0.325, 0.3], ["CLADDING", "HOLLOW", "FUEL"]):
    rect_cell.add(
        Region(Circle(radius=r), properties={PropertyType.MATERIAL: mat})
    )

# Display the cell's regions with the 'MATERIAL' property type colour map
rect_cell.show(PropertyType.MATERIAL)

# Build a lattice made of 36 Cartesian cells on 3 rings without central cell
lattice = CartesianLattice([], name='Cartesian Lattice')
lattice.add_rings_of_cells(rect_cell, 3)
# Apply a 1/8 symmetry
lattice.apply_symmetry(SymmetryType.EIGHTH)

# Display the lattice's regions with the 'MATERIAL' property type colour map
lattice.show(PropertyType.MATERIAL)

# Generate the output TDT file from the lattice
export_layout_to_tdt(
    lattice,
    os.path.join(
        os.path.dirname(sys.argv[0]), 'test_cartesian_assembly_eighth'
    ),
    TdtSetup(
        type_geo=LayoutGeometryType.RECTANGLE_EIGHT,
        symmetry_type=SymmetryType.EIGHTH
    )
)
