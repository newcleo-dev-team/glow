"""
Hexagonal assembly made of a central pin fuel cell sorrounded by 6 rings of
pin fuel cells. Symmetry is exploited by considering a one-sixth of the
complete hexagonal layout with rotational BCs on the interal borders.
"""
import os
import sys

from math import cos, pi

from glow.geometry_layouts.cells import *
from glow.geometry_layouts.lattices import *
from glow.main import TdtSetup, export_layout_to_tdt
from glow.generator.generator import *


# Declare the size of the edge for the hexagonal fuel cell
edge_length = 0.7852193995
# Build a hexagonal cell filled with 'COOLANT' and with four circular regions
fuel_cell = HexCell(
    side=edge_length,
    name="Hexagonal cell",
    base_props={PropertyType.MATERIAL: "COOLANT"}
)
for r, mat in zip(
        [0.525, 0.465, 0.45, 0.1],
        ["CLADDING", "HOLLOW", "FUEL", "HOLLOW"]
    ):
    fuel_cell.add(
        Region(Circle(radius=r), properties={PropertyType.MATERIAL: mat})
    )
# Rotate the cell by 90°
fuel_cell.rotate(90)
# Display the cell's regions with the 'MATERIAL' property type colour map
fuel_cell.show(PropertyType.MATERIAL)

# Build the lattice made of a central fuel cell surrounded by 6 rings of fuel
# cells around the central one
lattice = HexLattice([fuel_cell], name="Fuel lattice")
# Add a specific number of rings of cells to the lattice
lattice.add_rings_of_cells(fuel_cell, 6)

# Enclose the lattice in a hexagonal cell being its box
layer_thickness = 0.05
assembly = HexCell(
    side=lattice.dimensions[0] + 2*layer_thickness*cos(pi/3),
    base_props={PropertyType.MATERIAL: "COOLANT"},
    name="Fuel assembly"
)
assembly.add(
    Region(
        Hexagon(edge_length=lattice.dimensions[0] + layer_thickness*cos(pi/3))
        - lattice.shape,
        properties={PropertyType.MATERIAL: "CLADDING"}
    )
)
assembly.add(lattice)

# Apply the 'SIXTH' symmetry type
assembly.apply_symmetry(SymmetryType.SIXTH)
# Display the assembly's regions with the 'MATERIAL' property type colour map
assembly.show(PropertyType.MATERIAL)

# Generate the output TDT file from the layout
export_layout_to_tdt(
    assembly,
    os.path.join(
        os.path.dirname(sys.argv[0]), 'test_hex_assembly_sixth_rota'
    ),
    TdtSetup(
        type_geo=LayoutGeometryType.RA60,
        symmetry_type=SymmetryType.SIXTH
    )
)
